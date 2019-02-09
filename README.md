Venngraph Layout
================
d3.js plugin which renders venn diagram with edge information.
Based on [venn.js by benfred](https://github.com/benfred/venn.js).

## Motivation
This layout is basically same as venn diagram, but places node based on edge connection that makes two node closer if connected or in same group. Node with different group and no connection will be distantly placed.
Nelder-mead method is used to decide node position.

## Usage
```js
var nodes = [
    {group: 'A', name: 'n1'},       // n1 node is in group A
    {group: 'A__B', name: 'n2'},    // n2 node is in intersection of group A and B
    ...
];
var edges = [
    {source: 'n1', target: 'n2'},   // undirected edge information, which indicates connection between edge n1 and n2.
    ...
];
var vg = vennGraph(div).width(600).height(500).setNodes(nodes).setEdges(edges).prepare();
```

## Example
- [simple example](https://kuna.github.io/d3.venngraph.js/examples/example.html)

## License
MIT License