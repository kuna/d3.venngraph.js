/*
 * MIT License
 *
 * by @lazykuna
 * requires: d3.venn.js
 */


/*
 * @brief returns color-list given literal string.
 *        (ex: clr("A__B") == [clr("A"), clr("B")])
 * 
 * @param list    groups to be used.
 * @param preset  RGB color list to be used.
 */
var ColorMultiset = function(list, preset)
{
  var clrdict = {};
  var presetidx = 0;
  function color(s)
  {
    if (s.indexOf('__') > 0)
    {
      var a = s.split('__');
      var r = [];
      for (var i=0; i<a.length; i++)
      {
        r.push(color(a[i]));
      }
      return r;
    }

    if (!(s in clrdict))
    {
      clrdict[s] = preset[presetidx++];
      presetidx = presetidx % preset.length;
    }
    return clrdict[s];
  }

  // init
  if (!list) { list = []; }
  if (!preset) { preset = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]; }
  list.sort();
  for (var i=0; i<list.length; i++)
  {
    color(list[i]);
  }

  this.color = color;
  return this;
}


/*
 * @brief 
 * generates venn, node, edge objects based on node / edge data.
 * 
 * @param selection d3 selected object
 * @param clr       a dict which has group name as key and color string as value.
 *                  Automatically set if undefined.
 */
var vennGraph = function(selection, clr)
{
  var parameters = {
    width: 600,
    height: 500,
    maxIterations: 200,
    minErrorDelta: 1e-5,
    rho: 20,
    chi: 40,
    psi: -0.5,
    sigma: 0.5,
    maxNeighbors: 10
  };
  var vennChart = null;
  var vennData = null;
  var nodes = [];
  var edges = [];
  var nodes_pos = [];
  var edges_pos = [];
  var group_count = {};
  var circles = {};
  var exports = this;
  var svg = null;

  // nelderMead method
  // f: loss function
  // x0: initial value. re-calculates x0 to optimized value
  function weightedSum(ret, w1, v1, w2, v2) {
    for (var j = 0; j < ret.length; ++j) {
      ret[j] = w1 * v1[j] + w2 * v2[j];
    }
  }
  function nelderMead(f, x0, parameters) {
    parameters = parameters || {};
	
    var maxIterations = parameters.maxIterations || x0.length * 200,
        nonZeroDelta = parameters.nonZeroDelta || 1.05,
        zeroDelta = parameters.zeroDelta || 0.001,
        minErrorDelta = parameters.minErrorDelta || 1e-6,
        minTolerance = parameters.minErrorDelta || 1e-5,
        rho = (parameters.rho !== undefined) ? parameters.rho : 1,
        chi = (parameters.chi !== undefined) ? parameters.chi : 2,
        psi = (parameters.psi !== undefined) ? parameters.psi : -0.5,
        sigma = (parameters.sigma !== undefined) ? parameters.sigma : 0.5,
        maxDiff;

    // initialize simplex.
    var N = x0.length,
        simplex = new Array(N + 1);
    simplex[0] = x0;

    simplex[0].fx = f(x0);

    simplex[0].id = 0;
    for (var i = 0; i < N; ++i) {
        var point = x0.slice();
        point[i] = point[i] ? point[i] * nonZeroDelta : zeroDelta;
        simplex[i+1] = point;
        simplex[i+1].fx = f(point);
        simplex[i+1].id = i+1;
    }


    function updateSimplex(value) {
        for (var i = 0; i < value.length; i++) {
            simplex[N][i] = value[i];
        }
        simplex[N].fx = value.fx;
    }

    var sortOrder = function(a, b) { return a.fx - b.fx; };

    var centroid = x0.slice(),
        reflected = x0.slice(),
        contracted = x0.slice(),
        expanded = x0.slice();

    for (var iteration = 0; iteration < maxIterations; ++iteration) {
        simplex.sort(sortOrder);

        if (parameters.history) {
            // copy the simplex (since later iterations will mutate) and
            // sort it to have a consistent order between iterations
            var sortedSimplex = simplex.map(function (x) {
                var state = x.slice();
                state.fx = x.fx;
                state.id = x.id;
                return state;
            });
            sortedSimplex.sort(function(a,b) { return a.id - b.id; });

            parameters.history.push({x: simplex[0].slice(),
                                     fx: simplex[0].fx,
                                     simplex: sortedSimplex});
        }

        maxDiff = 0;
        for (i = 0; i < N; ++i) {
            maxDiff = Math.max(maxDiff, Math.abs(simplex[0][i] - simplex[1][i]));
        }

        if ((Math.abs(simplex[0].fx - simplex[N].fx) < minErrorDelta) &&
            (maxDiff < minTolerance)) {
            break;
        }

        // compute the centroid of all but the worst point in the simplex
        for (i = 0; i < N; ++i) {
            centroid[i] = 0;
            for (var j = 0; j < N; ++j) {
                centroid[i] += simplex[j][i];
            }
            centroid[i] /= N;
        }

        // reflect the worst point past the centroid  and compute loss at reflected
        // point
        var worst = simplex[N];
        weightedSum(reflected, 1+rho, centroid, -rho, worst);
        reflected.fx = f(reflected);

        // if the reflected point is the best seen, then possibly expand
        if (reflected.fx < simplex[0].fx) {
            weightedSum(expanded, 1+chi, centroid, -chi, worst);
            expanded.fx = f(expanded);
            if (expanded.fx < reflected.fx) {
                updateSimplex(expanded);
            }  else {
                updateSimplex(reflected);
            }
        }

        // if the reflected point is worse than the second worst, we need to
        // contract
        else if (reflected.fx >= simplex[N-1].fx) {
            var shouldReduce = false;

            if (reflected.fx > worst.fx) {
                // do an inside contraction
                weightedSum(contracted, 1+psi, centroid, -psi, worst);
                contracted.fx = f(contracted);
                if (contracted.fx < worst.fx) {
                    updateSimplex(contracted);
                } else {
                    shouldReduce = true;
                }
            } else {
                // do an outside contraction
                weightedSum(contracted, 1-psi * rho, centroid, psi*rho, worst);
                contracted.fx = f(contracted);
                if (contracted.fx < reflected.fx) {
                    updateSimplex(contracted);
                } else {
                    shouldReduce = true;
                }
            }

            if (shouldReduce) {
                // if we don't contract here, we're done
                if (sigma >= 1) break;

                // do a reduction
                for (i = 1; i < simplex.length; ++i) {
                    weightedSum(simplex[i], 1 - sigma, simplex[0], sigma, simplex[i]);
                    simplex[i].fx = f(simplex[i]);
                }
            }
        } else {
            updateSimplex(reflected);
        }
    }

    simplex.sort(sortOrder);
    return {fx : simplex[0].fx,
            x : simplex[0]};
  }

  // simple square root
  function loss_from_vector(node_vec, circles)
  {
    function sq(x) { return x*x; }
    function dist(a,b) { return Math.sqrt(sq(a.x-b.x)+sq(a.y-b.y)); }
    function dist_from_idx(a,b) { return Math.sqrt(sq(node_vec[a*2]-node_vec[b*2])+sq(node_vec[a*2+1]-node_vec[b*2+1])); }
    // closer to center of circles
    var dist_venn = 0;
    for (var i=0; i<nodes.length; i++)
    {
      dist_venn += dist(circles[nodes[i].group], {'x':node_vec[i*2], 'y':node_vec[i*2+1]});
    }
    // closer to neighbor (edge connected)
    var dist_neighbor = 0;
    for (var _=0; _<nodes.length; _++)
    {
      var np = nodes_pos[_];
      var i,x=0;
      for (i=0; i<np.neighbors.length && i<parameters.maxNeighbors; i++)
      {
        x += dist_from_idx(np.neighbors[i].idx, np.idx);
      }
      if (i > 0) x /= i;
      dist_neighbor += x;
    }
    // away from nearby nodes (in same group)
    var dist_nearby = 0;
    for (var _=0; _<nodes.length; _++)
    {
      var np = nodes_pos[_];
      var n = np.node;
      var i,x = 0;
      for (i=0; i<circles[n.group].nodes.length && i<parameters.maxNeighbors; i++)
      {
        var n2 = circles[n.group].nodes[i];
        x += dist_from_idx(n2.idx, n.idx);
      }
      if (i > 0) x /= i;
      dist_nearby += x;
    }

    var l = dist_nearby*0.5 - dist_neighbor - dist_venn*1.5;
    return -l;

  }

  function make_solution_vector(nodes)
  {
    var v = [];
    for (var i=0; i<nodes.length; i++)
    {
      v.push(nodes[i].x);
      v.push(nodes[i].y);
    }
    return v;
  }

  function fill_from_solution(nodes, solution)
  {
    for (var i=0; i<nodes.length; i++)
    {
      nodes[i].x = solution[2*i];
      nodes[i].y = solution[2*i+1];
    }
  }

  function prepare()
  {
    // prepare group count data
    group_count = {};
    for (var i=0; i<nodes.length; i++)
    {
      nodes[i].idx = i;   // additional information for vennGraph
      if (!(nodes[i].group in group_count)) group_count[nodes[i].group] = 0;
    }
    var kconds = []
    for (var k in group_count) { kconds.push( k.split('__') ); }
    for (var i=0; i<nodes.length; i++)
    {
      var cond = nodes[i].group;
      for (var j=0; j<kconds.length; j++)
      {
        var condvalid = true;
        for (var k=0; k<kconds[j].length; k++)
        {
          if (cond.indexOf(kconds[j][k])<0)
          {
            condvalid = false;
            break;
          }
        }
        if (condvalid) { group_count[Object.keys(group_count)[j]]++; }
      }
    }

    // prepare clr if undefined.
    if (clr === undefined)
    {
      clr = ColorMultiset(Object.keys(group_count)).color;
    }

    // initialize vennChart and fetch venn circle data
    var grouplist = [];
    for (var k in group_count)
    {
      grouplist.push( {'sets':k.split('__'), 'size':group_count[k], 'figure':group_count[k], 'label':k} );
    }
    vennChart = venn.VennDiagram().width(parameters.width).height(parameters.height).colours(clr);
    vennData = vennChart(selection.datum(grouplist));
    svg = selection.select('svg');  // fetch svg after vennChart
    setCircle(vennData.circles);

    console.log('calculating position of nodes. please wait ...');
    // initialize node position (random)
    for (var i=0; i<nodes.length; i++)
    {
      var ang = Math.random() * 2 * Math.PI;
      var r = circles[nodes[i].group].radius * Math.random() * 0.9;
      nodes[i].x = circles[nodes[i].group].x + r * Math.cos(ang);
      nodes[i].y = circles[nodes[i].group].y + r * Math.sin(ang);
    }

    // prepare neighbor information of node / edge
    nodes_pos = [];
    edges_pos = [];
    var nodedict = {};
    for (var i=0; i<nodes.length; i++)
    {
      var centre = circles[nodes[i].group];
      var ang_rand = Math.random() * 2 * Math.PI;
      var rad_rand = centre.radius * Math.random() * 0.9;
      var nm = {
        node: nodes[i],
        neighbors: [],
        x: centre.x + r * Math.cos(ang),
        y: centre.x + r * Math.cos(ang),
        idx: i
      };
      nodes_pos.push(nm);
      nodedict[nodes[i].name] = nm;
    }
    for (var i=0; i<edges.length; i++)
    {
      if (!(edges[i].source in nodedict)) throw("edge source '" + edges[i].source + "' is not in node!");
      if (!(edges[i].target in nodedict)) throw("edge target '" + edges[i].target + "' is not in node!");
      var cur_edge = {
        edge: edges[i],
        nsource: nodedict[edges[i].source],
        ntarget: nodedict[edges[i].target],
        x1:0, y1:0, x2:0, y2:0
      };
      cur_edge.nsource.neighbors.push(cur_edge.ntarget);
      cur_edge.ntarget.neighbors.push(cur_edge.nsource);
      edges_pos.push(cur_edge);
    }    

    // serialize node position to vector before nelderMead method
    var initial = make_solution_vector(nodes_pos);
    // Original nelderMead
    // - do gradient all of the node's position, based on nelderMead method
    var solution = nelderMead(
      function(values) {
        return loss_from_vector(values, circles);
     }, 
    initial, parameters);
    if (isNaN(solution.fx))
      console.warn("warning: nelderMead may failed to find proper position for nodes.");

    // update node / edge pos
    fill_from_solution(nodes_pos, solution.x);
    for (var i=0; i<edges_pos.length; i++)
    {
      edges_pos[i].x1 = edges_pos[i].nsource.x;
      edges_pos[i].y1 = edges_pos[i].nsource.y;
      edges_pos[i].x2 = edges_pos[i].ntarget.x;
      edges_pos[i].y2 = edges_pos[i].ntarget.y;
    }

    return exports;
  }
  
  function vennNode()
  {
    //fill_from_solution(nodes, initial);
    // place nodes
    var enter = svg.append('g')
      .attr('class', 'venn-nodes').selectAll('circle').data(nodes_pos)
      .enter()
      .append('circle')
      .attr('r', '2px')
      .attr('cx', function(d) { return d.x; })
      .attr('cy', function(d) { return d.y; })
      .style("fill", function (d) {
			  return clr(d.node.group);
		  });

    var entertext = svg.append('g')
      .attr('class', 'venn-nodes-text').selectAll('text').data(nodes_pos)
      .enter()
      .append('text')
      .attr('x', function(d) { return d.x+2; })
      .attr('y', function(d) { return d.y+3; })
      .text(function(d) { return d.node.name; });
    
    return {'nodes': nodes_pos,
      'enter': enter,
      'entertext': entertext
      };
  }

  function vennEdge()
  {
    // use nsource, ntarget attributes
    var enter = svg.append('g')
      .attr('class', 'venn-edges').selectAll('line').data(edges_pos)
      .enter()
      .append('line')
      .attr('x1', function (d) { return d.nsource.x; })
      .attr('y1', function (d) { return d.nsource.y; })
      .attr('x2', function (d) { return d.ntarget.x; })
      .attr('y2', function (d) { return d.ntarget.y; });
    return {
      'edges': edges_pos,
      'enter': enter
    };
  }

  function setWidth(v)
  {
    parameters.width = v;
    return exports;
  }

  function setHeight(v)
  {
    parameters.height = v;
    return exports;
  }

  function setNodes(o)
  {
    nodes = o;
    if (circles.length) makeCircleNode();
    return exports;
  }

  function setEdges(o)
  {
    edges = o;
    return exports;
  }

  function setCircle(o)
  {
    // do little process before setting circle
    circles = {};
    for (var k in o)
    {
      circles[k] = Object.assign({}, o[k], {'nodes':[]});
    }
    if (nodes.length) makeCircleNode();
    return exports;
  }

  function makeCircleNode()
  {
    for (var i=0; i<nodes.length; i++)
    {
      // make intersection circles (virtual circle)
      if (!(nodes[i].group in circles))
      {
        circles[nodes[i].group] = {'x':0, 'y':0, 'nodes':[]};
        var node_groups = nodes[i].group.split('__');
        var no_x=0, no_y=0, de=0;
        for (var j=0; j<node_groups.length; j++)
        {
          var r = circles[node_groups[j]].radius;
          no_x += circles[node_groups[j]].x*r;
          no_y += circles[node_groups[j]].y*r;
          de += r;
        }
        circles[nodes[i].group].x = no_x/de;
        circles[nodes[i].group].y = no_y/de;
        circles[nodes[i].group].radius = de/node_groups.length;
      }
      // fill circle.nodes
      circles[nodes[i].group].nodes.push(nodes[i]);
    }
  }

  function getChart()
  {
    return vennChart;
  }

  function getChartData()
  {
    return vennData;
  }

  this.vennNode = vennNode;
  this.vennEdge = vennEdge;
  this.getChart = getChart;
  this.getChartData = getChartData;
  this.width = setWidth;
  this.height = setHeight;
  this.setNodes = setNodes;
  this.setEdges = setEdges;
  //this.setCircle = setCircle; // -- changed to internal method as venn.js integrated.
  this.prepare = prepare;
  return this;
};


