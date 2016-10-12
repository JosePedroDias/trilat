var LM = require('ml-curve-fitting');
var Matrix = LM.Matrix;
var math = Matrix.algebra;



//Distance function. p is the guessed point.
var euclidean = function(t,p,c){
    var rows = t.rows;
    var result = new Matrix(t.rows, 1);
    for(var i=0;i<rows;i++){
       result[i][0] = Math.sqrt(Math.pow(t[i][0]-p[0][0],2)+Math.pow(t[i][1]-p[1][0],2));
    }

    return result;
};



// var data =  [[0.0, 0.0, 10.0], [10.0, 10.0, 10], [10.0, 0.0, 14.142135]];



/*
 * Data should be an array of 3 3-sized arrays with [x, y, dist]
 * Example follows...
 * Input:
 * [
 *  [ 0.0,  0.0, 10.0],
 *  [10.0, 10.0, 10.0],
 *  [10.0,  0.0, 14.142135]
 * ]
 *
 * Output:
 *
 */

function trilat(data, allowedDist) {
    var nbPoints = data.length;
    var t = math.matrix(nbPoints,2);//[1:Npnt]'; // independent variable
    var y_data = math.matrix(nbPoints, 1);

    for(var i=0;i<nbPoints;i++){
        t[i][0] = data[i][0];
        t[i][1] = data[i][1];
        y_data[i][0] = data[i][2];
    }

    var weight = [1];
    var opts = [ 2, 100, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1 ];
    var consts = [];
    
    var Xs = [ data[0][0], data[1][0], data[2][0] ];
    var Ys = [ data[0][1], data[1][1], data[2][1] ];
    var minX = Math.min.apply(Math, Xs);
    var minY = Math.min.apply(Math, Ys);
    var maxX = Math.max.apply(Math, Xs);
    var maxY = Math.max.apply(Math, Ys);
    var avgX = ( Xs[0] + Xs[1] + Xs[2] ) / 3;
    var avgY = ( Ys[0] + Ys[1] + Ys[2] ) / 3;
    var ad = allowedDist || 0;
    
    var p_init = math.matrix([[avgX], [avgY]]);
    var p_min = math.matrix([[minX-ad], [minY-ad]]);
    var p_max = math.matrix([[maxX+ad], [maxY-ad]]);
    
    // https://github.com/mljs/curve-fitting/blob/master/Documentation.md
    var p_fit = LM.optimize(euclidean,p_init,t,y_data,weight,-0.01,p_min,p_max,consts,opts);
    p_fit = p_fit.p;
    
    // euclidean(t,p_fit,consts)
    
    return [ p_fit[0][0], p_fit[1][0] ];
}

module.exports = trilat;
