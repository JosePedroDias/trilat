# trilat

Simple computation of trilateration using
matlab's Levenberg Marquardt curve-fitting algorithm.


## Common usage

```javascript
var trilat = require('trilat');

var input = [
//      X     Y     R
    [ 0.0,  0.0, 10.0],
    [10.0, 10.0, 10.0],
    [10.0,  0.0, 14.142135]
];

var output = trilat(input);
// [ 2.205170988086251e-7, 9.999999779478834 ]
```


## Usage without dependencies

If you want to use this in the browser without browserify/webpack/etc,
include the `dist.js` file, which exposes the symbol `trilat` as a global function.
