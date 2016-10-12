var trilat = require('./index');

var input = [
    [ 0.0,  0.0, 10.0],
    [10.0, 10.0, 10.0],
    [10.0,  0.0, 14.142135]
];

var output = trilat(input);

console.log(output);
