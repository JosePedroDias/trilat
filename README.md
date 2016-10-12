# trilat

Simple computation of trilateration using
matlab's Levenberg Marquardt curve-fitting algorithm.


## Usage

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
