//////////////////////////
//// CONSTANTS and general tools
///////////////////////////

import {type, isScalar, isArrayOfNumbers, error, laloparse,
        printNumber} from './lalolibbase.js';

export let LALOLIB_ERROR = "";
export const EPS = 2.2205e-16;

export function isZero(x) {
  return (Math.abs(x) < EPS);
}

export function isInteger(x) {
  return (Math.floor(x) === x);
}

export let TICTOCstartTime;

export function tic(T) {
  if (typeof (TICTOCstartTime) === "undefined")
    TICTOCstartTime = [];
  if (typeof (T) === "undefined")
    T = 0;
  TICTOCstartTime[T] = new Date();
}

export function toc(T) {
  if (typeof (T) === "undefined") {
    T = 0;
  }

  if (typeof (TICTOCstartTime) !== "undefined" && typeof (TICTOCstartTime[T]) !== "undefined") {
    // Computing time
    let startTime = TICTOCstartTime[T];
    let endTime = new Date();
    let time = (endTime - startTime)/1000;  // in seconds
    return time;
  } else
    return undefined;
}


/**
 * @param {Float64Array}
 * @return {string}
 */
export function printVector(x) {
  const n = x.length;
  let str = "[ ";
  let i = 0;
  while (i < n - 1 && i < 5) {
    str += (isInteger(x[i]) ? x[i] : x[i].toFixed(3)) + "; ";
    i++;
  }
  if (i === n - 1)
    str += (isInteger(x[i]) ? x[i] : x[i].toFixed(3)) + " ]";
  else
    str += "... ] (length = " + n + ")";

  return str;
}


//////////////////////////////
// Matrix/vector creation
//////////////////////////////
/**
 * @constructor
 * @struct
 */

export class Matrix {
  constructor(m, n, values) {

    /** @const */ this.length = m;
    /** @const */ this.m = m;
    /** @const */ this.n = n;
    /** @const */ this.size = [m, n];
    /** @const */ this.type = "matrix";

    if (arguments.length === 2)
      this.val = new Float64Array(m*n); // simple m x n zeros
    else if (arguments.length === 3)
      this.val = new Float64Array(values); // m x n filled with values with copy
    else if (arguments.length === 4)
      this.val = values; // m x n filled with values without copy
  }

  get(i, j) {
    return this.val[i*this.n + j];
  }

  set(i, j, v) {
    this.val[i*this.n + j] = v;
  }

  /**
   * return a pointer-like object on a row in a matrix, not a copy!
   * @param {number}
   * @return {Float64Array}
   */
  row(i) {
    return this.val.subarray(i*this.n, (i + 1)*this.n);
  }

  /**
   * return a copy of the matrix as an Array of Arrays
   * (do not do this with too many rows...)
   * @return {Array}
   */
  toArray() {
    let A = new Array(this.m);
    let ri = 0;
    for (let i = 0; i < this.m; i++) {
      A[i] = new Array(this.n);
      for (let j = 0; j < this.n; j++) {
        A[i][j] = this.val[ri + j];
      }
      ri += this.n;
    }
    return A;
  }

  /**
   * return a view (not a copy) on the matrix as an Array of Float64Array
   * (do not do this with too many rows...)
   * @return {Array}
   */
  toArrayOfFloat64Array() {
    let A = new Array(this.m);
    for (let i = 0; i < this.m; i++) {
      A[i] = this.val.subarray(i*this.n, (i + 1)*this.n);
    }

    return A;
  }
}

export function array2mat(A) {
  return mat(A, true);
}

export function array2vec(a) {
  return vectorCopy(a);
}

export function vec2array(a) {
  return Array.apply([], a);
}

export function size(A, sizealongdimension) {
  let s;
  switch (type(A)) {
    case "string":
    case "boolean":
    case "number":
    case "Complex":
      s = [1, 1];
      break;
    case "vector":
    case "spvector":
    case "ComplexVector":
      s = [A.length, 1];
      break;
    case "matrix":
    case "spmatrix":
    case "ComplexMatrix":
      s = A.size;
      break;
    case "object":
      s = [1, 1];
      break;
    default:
      s = [1, 1];
      //error( "Cannot determine size of object" );
      break;
  }

  if (typeof (sizealongdimension) === "undefined")
    return s;
  else
    return s[sizealongdimension - 1];

}

export function ones(rows, cols) {
  // Create a matrix or vector full of ONES
  if (arguments.length === 1 || cols === 1) {
    let v = new Float64Array(rows);
    for (let i = 0; i < rows; i++) {
      v[i] = 1;
    }
    return v;
  } else {
    let M = new Matrix(rows, cols);
    const mn = rows*cols;
    for (let i = 0; i < mn; i++) {
      M.val[i] = 1;
    }
    return M;
  }
}

// Use zeros( m, n)
export function zeros(rows, cols) {
  // Create a matrix or vector of ZERO
  if (arguments.length === 1 || cols === 1) {
    return new Float64Array(rows);
  } else {
    return new Matrix(rows, cols);
  }
}

export function eye(m, n=m) {
  if (m === 1 && n === 1)
    return 1;

  let I = zeros(m, n);
  const e = (m < n) ? m : n;
  for (let i = 0; i < e; i++) {
    I.val[i*(n + 1)] = 1;
  }

  return I;
}

export function diag(A) {
  let i;
  let typeA = type(A);
  if (typeA === "vector") {
    let M = zeros(A.length, A.length);
    let j = 0;
    const stride = A.length + 1;
    for (i = 0; i < A.length; i++) {
      M.val[j] = A[i];
      j += stride;
    }
    return M;
  } else if (typeA === "matrix") {
    let n = Math.min(A.m, A.n);
    let v = new Float64Array(n);
    let j = 0;
    const stride2 = A.n + 1;
    for (i = 0; i < n; i++) {
      v[i] = A.val[j];
      j += stride2;
    }
    return v;
  } else if (typeA === "ComplexVector") {
    let M = new ComplexMatrix(A.length, A.length);
    let j = 0;
    const stride = A.length + 1;
    for (i = 0; i < A.length; i++) {
      M.re[j] = A.re[i];
      M.im[j] = A.im[i];
      j += stride;
    }
    return M;
  } else if (typeA === "ComplexMatrix") {
    let n = Math.min(A.m, A.n);
    let v = new ComplexVector(n);
    let j = 0;
    const stride2 = A.n + 1;
    for (i = 0; i < n; i++) {
      v.re[i] = A.re[j];
      v.im[i] = A.im[j];
      j += stride2;
    }
    return v;
  }
}

/**
 * @param {Matrix}
 * @return {Float64Array}
 */
export function vec(A) {
  return new Float64Array(A.val);
}

export function matrixCopy(A) {
  let t = type(A);
  switch (t) {
    case "vector":
      return vectorCopy(A);
    case "ComplexVector":
      return new ComplexVector(A);
    case "matrix":
      return new Matrix(A.m, A.n, A.val);
    case "ComplexMatrix":
      return new ComplexMatrix(A);
    case "Array":
      return arrayCopy(A);
    case "spvector":
    case "spmatrix":
      return A.copy();
    default:
      error("Error in matrixCopy(A): A is not a matrix nor a vector.");
      return undefined;
  }
}

/**
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function vectorCopy(a) {
  return new Float64Array(a);
}

/** Vector copy into another existing vector ( y = x )
 * (saves memory allocation)
 * @param {Float64Array}
 * @param {Float64Array}
 */
export function vectorCopyInto(x, y) {
  y.set(x);
}

/**
 * @param {Array}
 * @return {Array}
 */
export function arrayCopy(A) {
  let res = new Array(A.length);
  for (let i = 0; i < A.length; i++) {
    if (isScalar(A[i]))
      res[i] = A[i];	//does not copy 2D Arrays...
    else
      res[i] = matrixCopy(A[i]);
  }
  return res;
}

/**
 * Return enlarged matrix with one more row of zeros
 * NOTE: both matrices share the same storage and should not be used independently
 * so better use: A = appendRow(A); or just appendRow(A);
 * @param{Matrix}
 * @return{Matrix}
 */
export function appendRow(A) {
  let Aa = zeros(A.m + 1, A.n);
  Aa.val.set(A.val);
  return Aa;
}

/**
 * Reshape the dimensions of a vector/matrix
 * @param{{Float64Array|Matrix}}
 * @param{number}
 * @param{number}
 * @return{{Float64Array|Matrix}}
 */
export function reshape(A, m, n) {
  let R = undefined;
  let tA = type(A);
  if (tA === "vector") {
    if (m*n !== A.length) {
      error("Error in reshape(a,m,n): a.length = " + A.length + " !== m*n");
    } else {
      R = new Matrix(m, n, A);
    }
  } else if (tA === "matrix") {
    if (m*n !== A.m*A.n) {
      error("Error in reshape(A,m,n): A.m * A.n = " + A.m*A.n + " !== m*n");
    } else {
      if (n === 1)
        R = vectorCopy(A.val);
      else
        R = new Matrix(m, n, A.val);
    }
  } else
    error("Error in reshape(A): A is neither a vector nor a matrix.");
  return R;
}


////////////////////////
// slicing functions
////////////////////////

/*
	GET function : returns a copy of a subset of entries

	For MATRICES:

	get ( M, rows, cols ) => submatrix of M
	get ( M, rows ) 	  => subset of rows from M (equiv to rows(M,rows) )
	get ( M, [], cols )   => subset of cols (equiv to cols(M, cols) )
	get ( M, i, j)		  => M[i][j] converted to dense format (0 instead of undefined)
	get ( M ) 			  => M in dense format  (with 0 instead of undefined)

	For VECTORS:

	get ( v, rows ) 	  => subvector from v (equiv to rows(v,rows) )
	get ( v, i )		  => v[i] converted to dense format (0 instead of undefined)
	get ( v ) 			  => v in dense format  (with 0 instead of undefined)

*/
export function get(A, rowsrange, colsrange) {

  let typerows = typeof (rowsrange);
  let typecols = typeof (colsrange);

  if (arguments.length === 1)
    return matrixCopy(A);

  let typeA = type(A);
  if (typeA === "vector") {

    if (typerows === "number") {
      if (rowsrange >= 0 && rowsrange < A.length)
        return A[rowsrange];	// get v[i]
      else {
        error("Error in a[i] = get(a,i): Index i=" + rowsrange + " out of bounds [0," + (A.length - 1) + "]");
        return undefined;
      }
    } else {
      return getSubVector(A, rowsrange);
    }
  } else if (typeA === "matrix") {

    if (typerows === "number")
      rowsrange = [rowsrange];

    if (typecols === "number")
      colsrange = [colsrange];

    if (rowsrange.length === 1 && colsrange.length === 1)
      return A.val[rowsrange[0]*A.n + colsrange[0]];	// get ( A, i, j)

    if (rowsrange.length === 0)
      return getCols(A, colsrange);// get(A,[],4) <=> cols(A,4)

    if (colsrange.length === 0)
      return getRows(A, rowsrange);// get(A,3,[]) <=> rows(A,3)

    // otherwise:
    return getSubMatrix(A, rowsrange, colsrange);

  } else if (typeA === "Array") {
    if (typerows === "number")
      return A[rowsrange];
    else
      return getSubArray(A, rowsrange);
  } else if (typeA === "spmatrix") {

    if (typerows === "number")
      rowsrange = [rowsrange];

    if (typecols === "number")
      colsrange = [colsrange];

    if (rowsrange.length === 1 && colsrange.length === 1)
      return A.get(rowsrange[0], colsrange[0]);   // get ( A, i, j)

    if (rowsrange.length === 1 && A.rowmajor)
      return A.row(rowsrange[0]);
    if (colsrange.length === 1 && !A.rowmajor)
      return A.col(colsrange[0]);

    if (colsrange.length === 0)
      return spgetRows(A, rowsrange);
    if (rowsrange.length === 0)
      return spgetCols(A, colsrange);

    // TODO
  } else if (typeA === "spvector") {

    if (typerows === "number")
      return A.get(rowsrange);	// get v[i]
    else
      return getSubspVector(A, rowsrange);//TODO
  } else if (typeA === "ComplexVector") {
    if (typerows === "number")
      return A.get(rowsrange);	// get v[i]
    else
      return A.getSubVector(rowsrange);
  } else if (typeA === "ComplexMatrix") {

    if (typerows === "number")
      rowsrange = [rowsrange];

    if (typecols === "number")
      colsrange = [colsrange];

    if (rowsrange.length === 1 && colsrange.length === 1)
      return A.get(i, j);

    if (rowsrange.length === 0)
      return A.getCols(colsrange);// get(A,[],4) <=> cols(A,4)

    if (colsrange.length === 0)
      return A.getRows(rowsrange);// get(A,3,[]) <=> rows(A,3)

    // otherwise:
    return A.getSubMatrix(rowsrange, colsrange);
  }
  return undefined;
}

export function getSubMatrix(A, rowsrange, colsrange) {
  let n = colsrange.length;
  let i;
  let j;
  let res;
  if (n === 1) {
    res = new Float64Array(rowsrange.length);
    for (i = 0; i < rowsrange.length; i++) {
      res[i] = A.val[rowsrange[i]*A.n + colsrange[0]];
    }
  } else {
    res = new Matrix(rowsrange.length, n);
    let r = 0;

    for (i = 0; i < rowsrange.length; i++) {
      let rA = rowsrange[i]*A.n;
      for (j = 0; j < n; j++) {
        res.val[r + j] = A.val[rA + colsrange[j]];
      }
      r += n;
    }
  }
  return res;
}

export function getRows(A, rowsrange) {
  let n = rowsrange.length;
  if (n > 1) {
    let res = new Matrix(n, A.n);
    let r = 0;
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < A.n; j++) {
        res.val[r + j] = A.val[rowsrange[i]*A.n + j];
      }
      r += A.n;
    }
    return res;
  } else
    return vectorCopy(A.val.subarray(rowsrange[0]*A.n, rowsrange[0]*A.n + A.n));
}

export function getCols(A, colsrange) {
  let m = A.m;
  let n = colsrange.length;
  if (n > 1) {
    let res = new Matrix(m, n);
    let r = 0;
    let rA = 0;
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < n; j++) {
        res.val[r + j] = A.val[rA + colsrange[j]];
      }

      r += n;
      rA += A.n;
    }
    return res;
  } else {
    let res = new Float64Array(m);
    let r = 0;
    for (let i = 0; i < m; i++) {
      res[i] = A.val[r + colsrange[0]];
      r += A.n;
    }
    return res;
  }
}

/**
 * @param {Float64Array}
 * @param {Array}
 * @return {Float64Array}
 */
export function getSubVector(a, rowsrange) {
  const n = rowsrange.length;
  let res = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    res[i] = a[rowsrange[i]];
  }
  return res;
}

/**
 * @param {Array}
 * @param {Array}
 * @return {Array}
 */
export function getSubArray(a, rowsrange) {
  const n = rowsrange.length;
  let res = new Array(n);
  for (let i = 0; i < n; i++) {
    res[i] = a[rowsrange[i]];
  }
  return res;
}


export function getrowref(A, i) {
  // return a pointer-like object on a row in a matrix, not a copy!
  return A.val.subarray(i*A.n, (i + 1)*A.n);
}

/*
	SET function : set values in a subset of entries of a matrix or vector

	For MATRICES:

	set ( M, rows, cols, A ) => submatrix of M = A
	set ( M, rows, A ) 	     => subset of rows from M = A
	set ( M, [], cols, A )   => subset of cols from M = A
	set ( M, i, [], A )   	 => fill row M[i] with vector A (transposed)
	set ( M, i, j, A)	     => M[i][j] = A

	For VECTORS:

	set ( v, rows, a ) 	  => subvector from v = a
	set ( v, i , a)		  => v[i] = a

*/
export function set(A, rowsrange, colsrange, B) {
  let i;
  let j;
  let k;
  let l;
  let n;

  let typerows = typeof (rowsrange);
  let typecols = typeof (colsrange);

  if (arguments.length === 1)
    return undefined;

  let typeA = type(A);
  if (typeA === "vector") {
    B = colsrange;
    if (typerows === "number") {
      A[rowsrange] = B;
      return B;
    } else if (rowsrange.length === 0)
      rowsrange = range(A.length);

    if (size(B, 1) === 1) {
      setVectorScalar(A, rowsrange, B);
    } else {
      setVectorVector(A, rowsrange, B);
    }
    return B;
  } else if (typeA === "matrix") {

    if (typerows === "number")
      rowsrange = [rowsrange];
    if (typecols === "number")
      colsrange = [colsrange];

    if (rowsrange.length === 1 && colsrange.length === 1) {
      A.val[rowsrange[0]*A.n + colsrange[0]] = B;
      return B;
    }

    if (rowsrange.length === 0) {
      setCols(A, colsrange, B);
      return B;
    }

    if (colsrange.length === 0) {
      setRows(A, rowsrange, B);
      return B;
    }

    // Set a submatrix
    let sB = size(B);
    let tB = type(B);
    if (sB[0] === 1 && sB[1] === 1) {
      if (tB === "number")
        setMatrixScalar(A, rowsrange, colsrange, B);
      else if (tB === "vector")
        setMatrixScalar(A, rowsrange, colsrange, B[0]);
      else
        setMatrixScalar(A, rowsrange, colsrange, B.val[0]);
    } else {
      if (colsrange.length === 1)
        setMatrixColVector(A, rowsrange, colsrange[0], B);
      else if (rowsrange.length === 1) {
        if (tB === "vector")
          setMatrixRowVector(A, rowsrange[0], colsrange, B);
        else
          setMatrixRowVector(A, rowsrange[0], colsrange, B.val);
      } else
        setMatrixMatrix(A, rowsrange, colsrange, B);
    }
    return B;
  } else if (typeA === "ComplexVector") {
    B = colsrange;
    if (typerows === "number") {
      A.set(rowsrange, B);
      return B;
    } else if (rowsrange.length === 0)
      rowsrange = range(A.length);

    if (size(B, 1) === 1) {
      A.setVectorScalar(rowsrange, B);
    } else {
      A.setVectorVector(rowsrange, B);
    }
    return B;
  }
}

export function setVectorScalar(A, rowsrange, B) {
  let i;
  for (i = 0; i < rowsrange.length; i++) {
    A[rowsrange[i]] = B;
  }
}

export function setVectorVector(A, rowsrange, B) {
  let i;
  for (i = 0; i < rowsrange.length; i++) {
    A[rowsrange[i]] = B[i];
  }
}

export function setMatrixScalar(A, rowsrange, colsrange, B) {
  let i;
  let j;
  let m = rowsrange.length;
  let n = colsrange.length;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      A.val[rowsrange[i]*A.n + colsrange[j]] = B;
    }
  }
}

export function setMatrixMatrix(A, rowsrange, colsrange, B) {
  let i;
  let j;
  let m = rowsrange.length;
  let n = colsrange.length;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      A.val[rowsrange[i]*A.n + colsrange[j]] = B.val[i*B.n + j];
    }
  }
}

export function setMatrixColVector(A, rowsrange, col, B) {
  let i;
  let m = rowsrange.length;
  for (i = 0; i < m; i++) {
    A.val[rowsrange[i]*A.n + col] = B[i];
  }
}

export function setMatrixRowVector(A, row, colsrange, B) {
  let j;
  let n = colsrange.length;
  for (j = 0; j < n; j++) {
    A.val[row*A.n + colsrange[j]] = B[j];
  }
}

export function setRows(A, rowsrange, B) {
  let i;
  let j;
  let m = rowsrange.length;
  let rA;
  switch (type(B)) {
    case "vector":
      for (i = 0; i < m; i++) {
        rA = rowsrange[i]*A.n;
        for (j = 0; j < B.length; j++) {
          A.val[rA + j] = B[j];
        }
      }
      break;
    case "matrix":
      let rB = 0;
      for (i = 0; i < m; i++) {
        rA = rowsrange[i]*A.n;
        for (j = 0; j < B.n; j++) {
          A.val[rA + j] = B.val[rB + j];
        }
        rB += B.n;
      }
      break;
    default:
      for (i = 0; i < m; i++) {
        rA = rowsrange[i]*A.n;
        for (j = 0; j < A.n; j++) {
          A.val[rA + j] = B;
        }
      }
      break;
  }
}

export function setCols(A, colsrange, B) {
  let i;
  let m = A.m;
  let n = colsrange.length;
  let r = 0;
  switch (type(B)) {
    case "vector":
      for (i = 0; i < m; i++) {
        for (let j = 0; j < n; j++) {
          A.val[r + colsrange[j]] = B[i];
        }
        r += A.n;
      }
      break;
    case "matrix":
      for (i = 0; i < m; i++) {
        for (let j = 0; j < n; j++) {
          A.val[r + colsrange[j]] = B.val[i*B.n + j];
        }
        r += A.n;
      }
      break;
    default:
      for (i = 0; i < m; i++) {
        for (let j = 0; j < n; j++) {
          A.val[r + colsrange[j]] = B;
        }
        r += A.n;
      }
      break;
  }
}

export function dense(A) {
  return A;
}

// Support
export function supp(x) {
  const tx = type(x);
  if (tx === "vector") {
    let indexes = [];
    let i;
    for (i = 0; i < x.length; i++) {
      if (!isZero(x[i]))
        indexes.push(i);
    }

    return indexes;
  } else if (tx === "spvector") {
    return new Float64Array(x.ind);
  } else
    return undefined;
}

// Range
export function range(start, end, inc=1) {
  // python-like range function
  // returns [0,... , end-1]
  if (typeof (start) === "undefined")
    return [];

  if (typeof (end) === "undefined") {
    let end = start;
    start = 0;
  }

  if (start === end - inc) {
    return start;
  } else if (start === end) {
    return [];
  } else if (start > end) {
    if (inc > 0)
      inc *= -1;
    let r = new Array(Math.floor((start - end)/Math.abs(inc)));
    let k = 0;
    for (let i = start; i > end; i += inc) {
      r[k] = i;
      k++;
    }
  } else {
    let r = new Array(Math.floor((end - start)/inc));
    let k = 0;
    for (let i = start; i < end; i += inc) {
      r[k] = i;
      k++;
    }
  }
  return r;
}

// Swaping
/**
 * @param {Matrix}
 */
export function swaprows(A, i, j) {
  if (i !== j) {
    let ri = i*A.n;
    let rj = j*A.n;
    let tmp = vectorCopy(A.val.subarray(ri, ri + A.n));
    A.val.set(vectorCopy(A.val.subarray(rj, rj + A.n)), ri);
    A.val.set(tmp, rj);
  }
}

/**
 * @param {Matrix}
 */
export function swapcols(A, j, k) {
  if (j !== k) {
    let tmp = getCols(A, [j]);
    setCols(A, [j], getCols(A, [k]));
    setCols(A, [k], tmp);
  }
}

//////////////////////////
// Random numbers
////////////////////////////

// Gaussian random number (mean = 0, variance = 1;
//	Gaussian noise with the polar form of the Box-Muller transformation
export function randnScalar() {

  let x1;
  let x2;
  let w;
  let y1;
  let y2;
  do {
    x1 = 2.0*Math.random() - 1.0;
    x2 = 2.0*Math.random() - 1.0;
    w = x1*x1 + x2*x2;
  } while (w >= 1.0);

  w = Math.sqrt((-2.0*Math.log(w))/w);
  y1 = x1*w;
  y2 = x2*w;

  return y1;
}

export function randn(dim1, dim2) {
  let res;

  if (typeof (dim1) === "undefined" || (dim1 === 1 && typeof (dim2) == "undefined") || (dim1 === 1 && dim2 == 1)) {
    return randnScalar();
  } else if (typeof (dim2) === "undefined" || dim2 === 1) {
    res = new Float64Array(dim1);
    for (let i = 0; i < dim1; i++) {
      res[i] = randnScalar();
    }

    return res;
  } else {
    res = zeros(dim1, dim2);
    for (let i = 0; i < dim1*dim2; i++) {
      res.val[i] = randnScalar();
    }
    return res;
  }
}

// Uniform random numbers
/*
 * @param{number}
 * @return{Float64Array}
 */
export function randVector(dim1) {
  let res = new Float64Array(dim1);
  for (let i = 0; i < dim1; i++) {
    res[i] = Math.random();
  }
  return res;
}

/*
 * @param{number}
 * @param{number}
 * @return{Matrix}
 */
export function randMatrix(dim1, dim2) {
  const n = dim1*dim2;
  let res = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    res[i] = Math.random();
  }
  return new Matrix(dim1, dim2, res, true);
}

export function rand(dim1, dim2) {
  let res;
  if (typeof (dim1) === "undefined" || (dim1 === 1 && typeof (dim2) == "undefined") || (dim1 === 1 && dim2 == 1)) {
    return Math.random();
  } else if (typeof (dim2) === "undefined" || dim2 === 1) {
    return randVector(dim1);
  } else {
    return randMatrix(dim1, dim2);
  }
}

export function randnsparse(NZratio, dim1, dim2) {
  // Generates a sparse random matrix with NZratio * dim1*dim2 (or NZ if NZratio > 1 ) nonzeros
  let NZ;
  if (NZratio > 1)
    NZ = NZratio;
  else
    NZ = Math.floor(NZratio*dim1*dim2);

  let indexes;
  let i;
  let j;
  let k;
  let res;

  if (typeof (dim1) === "undefined") {
    return randn();
  } else if (typeof (dim2) === "undefined" || dim2 === 1) {

    indexes = randperm(dim1);

    res = zeros(dim1);
    for (i = 0; i < NZ; i++) {
      res[indexes[i]] = randn();
    }
    return res;
  } else {
    res = zeros(dim1, dim2);
    indexes = randperm(dim1*dim2);
    for (k = 0; k < NZ; k++) {
      i = Math.floor(indexes[k]/dim2);
      j = indexes[k] - i*dim2;
      res.val[i*dim2 + j] = randn();
    }
    return res;
  }
}

export function randsparse(NZratio, dim1, dim2) {
  // Generates a sparse random matrix with NZratio * dim1*dim2 (or NZ if NZratio > 1 ) nonzeros
  if (typeof (dim2) === "undefined")
    dim2 = 1;

  let NZ;
  if (NZratio > 1)
    NZ = NZratio;
  else
    NZ = Math.floor(NZratio*dim1*dim2);

  let indexes;
  let i;
  let j;
  let k;
  let res;

  if (typeof (dim1) === "undefined") {
    return randn();
  } else if (dim2 === 1) {

    indexes = randperm(dim1);

    res = zeros(dim1);
    for (i = 0; i < NZ; i++) {
      res[indexes[i]] = Math.random();
    }
    return res;
  } else {
    res = zeros(dim1, dim2);
    indexes = randperm(dim1*dim2);

    for (k = 0; k < NZ; k++) {
      i = Math.floor(indexes[k]/dim2);
      j = indexes[k] - i*dim2;
      res.val[i*dim2 + j] = Math.random();
    }
    return res;
  }
}

export function randperm(x) {
  // return a random permutation of x (or of range(x) if x is a number)

  if (typeof (x) === "number") {
    let perm = range(x);
  } else {
    let perm = new Float64Array(x);
  }
  let i;
  let j;
  let k;

  // shuffle
  for (i = perm.length - 1; i > 1; i--) {
    j = Math.floor(Math.random()*i);
    k = perm[j];
    perm[j] = perm[i];
    perm[i] = k;
  }
  return perm;
}

///////////////////////////////
/// Basic Math function: give access to Math.* JS functions
///  and vectorize them
///////////////////////////////

import {spMatrix} from './sparse.js';

// automatically generate (vectorized) wrappers for Math functions
export let MathFunctions = Object.getOwnPropertyNames(Math);
export const VecMath = {};

for (let mf in MathFunctions) {
  let globalPref = "VecMath.";
  let key = MathFunctions[mf];

  if (typeof Math[key] === "function") {
    if (eval("Math." + MathFunctions[mf] + ".length") === 1) {
      // this is a function of a scalar
      // make generic function:
      eval(globalPref + MathFunctions[mf] + " = function (x) { return apply(Math." + MathFunctions[mf] + " , x );};");
      // make vectorized version:
      eval(globalPref + MathFunctions[mf] + "Vector = function (x) { return applyVector(Math." + MathFunctions[mf] + " , x );};");
      // make matrixized version:
      eval(globalPref + MathFunctions[mf] + "Matrix = function (x) { return applyMatrix(Math." + MathFunctions[mf] + " , x );};");
    }
  } else if (typeof Math[key] === "number") {
    // Math constant:
    eval(globalPref + MathFunctions[mf] + " = Math." + MathFunctions[mf]);
  }
}

/*
let s = 'let {';
for (let k in VecMath) {
  s += `${k}, `;
}
s += '} = VecMath';

console.log(s);
*/

let {abs, absVector, absMatrix, acos, acosVector, acosMatrix, acosh, acoshVector
      , acoshMatrix, asin, asinVector, asinMatrix, asinh, asinhVector, asinhMatrix,
      atan, atanVector, atanMatrix, atanh, atanhVector, atanhMatrix, ceil, ceilVector,
      ceilMatrix, cbrt, cbrtVector, cbrtMatrix, expm1, expm1Vector, expm1Matrix, clz32,
      clz32Vector, clz32Matrix, cos, cosVector, cosMatrix, cosh, coshVector, coshMatrix,
      exp, expVector, expMatrix, floor, floorVector, floorMatrix, fround, froundVector,
      froundMatrix, log, logVector, logMatrix, log1p, log1pVector, log1pMatrix,
      log2, log2Vector, log2Matrix, log10, log10Vector, log10Matrix, round, roundVector,
      roundMatrix, sign, signVector, signMatrix, sin, sinVector, sinMatrix, sinh, sinhVector,
      sinhMatrix, sqrt, sqrtVector, sqrtMatrix, tan, tanVector, tanMatrix, tanh, tanhVector,
      tanhMatrix, trunc, truncVector, truncMatrix, E, LN10, LN2, LOG10E, LOG2E, PI, SQRT1_2,
      SQRT2, fract, fractVector, fractMatrix, tent, tentVector,
      tentMatrix, } = VecMath;

export function apply(f, x) {
  // Generic wrapper to apply scalar functions
  // element-wise to vectors and matrices
  if (typeof (f) !== "function")
    return undefined;
  switch (type(x)) {
    case "number":
      return f(x);
      break;
    case "Complex":
      let ComplexFunctions = ["exp", "abs"];
      let fc = ComplexFunctions.indexOf(f.name);
      if (fc >= 0)
        return eval(ComplexFunctions[fc] + "Complex(x);");
      else {
        error("This function has no Complex counterpart (yet).");
        return undefined;
      }
      break;
    case "vector":
      return applyVector(f, x);
      break;
    case "spvector":
      return applyspVector(f, x);
      break;
    case "ComplexVector":
      if (f.name === "abs")
        return absComplex(x);
      else
        return applyComplexVector(f, x);
      break;
    case "matrix":
      return applyMatrix(f, x);
      break;
    case "spmatrix":
      return applyspMatrix(f, x);
      break;
    case "ComplexMatrix":
      if (f.name === "abs")
        return absComplex(x);
      else
        return applyComplexMatrix(f, x);
      break;
    default:
      return "undefined";
  }
}

export function applyVector(f, x) {
  const nv = x.length;
  let res = new Float64Array(nv);
  for (let i = 0; i < nv; i++) {
    res[i] = f(x[i]);
  }
  return res;
}

export function applyComplexVector(f, x) {
  const nv = x.length;
  let res = new ComplexVector(nv);
  for (let i = 0; i < nv; i++) {
    res.set(i, f(x.get(i)));
  }
  return res;
}

export function applyComplexMatrix(f, x) {
  const m = x.m;
  const n = x.n;
  let res = new ComplexMatrix(m, n);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      res.set(i, j, f(x.get(i, j)));
    }
  }
  return res;
}

export function applyMatrix(f, x) {
  return new Matrix(x.m, x.n, applyVector(f, x.val), true);
}

///////////////////////////////
/// Operators
///////////////////////////////

export function mul(a, b) {
  let sa = size(a);
  let sb = size(b);
  if (!isScalar(a) && sa[0] === 1 && sa[1] === 1)
    a = get(a, 0, 0);
  if (!isScalar(b) && sb[0] === 1 && sb[1] === 1)
    b = get(b, 0, 0);

  switch (type(a)) {
    case "number":
      switch (type(b)) {
        case "number":
          return a*b;
          break;
        case "Complex":
          return mulComplexReal(b, a);
          break;
        case "vector":
          return mulScalarVector(a, b);
          break;
        case "spvector":
          return mulScalarspVector(a, b);
          break;
        case "ComplexVector":
          return mulScalarComplexVector(a, b);
          break;
        case "matrix":
          return mulScalarMatrix(a, b);
          break;
        case "spmatrix":
          return mulScalarspMatrix(a, b);
          break;
        case "ComplexMatrix":
          return mulScalarComplexMatrix(a, b);
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "Complex":
      switch (type(b)) {
        case "number":
          return mulComplexReal(a, b);
          break;
        case "Complex":
          return mulComplex(a, b);
          break;
        case "vector":
          return mulComplexVector(a, b);
          break;
        case "ComplexVector":
          return mulComplexComplexVector(a, b);
          break;
        case "spvector":
          return mulComplexspVector(a, b);
          break;
        case "matrix":
          return mulComplexMatrix(a, b);
          break;
        case "ComplexMatrix":
          return mulComplexComplexMatrix(a, b);
          break;
        case "spmatrix":
          return mulComplexspMatrix(a, b);
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "vector":
      switch (type(b)) {
        case "number":
          return mulScalarVector(b, a);
          break;
        case "Complex":
          return mulComplexVector(b, a);
          break;
        case "vector":
          if (a.length !== b.length) {
            error("Error in mul(a,b) (dot product): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return dot(a, b);
          break;
        case "spvector":
          if (a.length !== b.length) {
            error("Error in mul(a,b) (dot product): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return dotspVectorVector(b, a);
          break;
        case "ComplexVector":
          if (a.length !== b.length) {
            error("Error in mul(a,b) (dot product): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return dotComplexVectorVector(b, a);
          break;
        case "matrix":
          if (b.m === 1)
            return outerprodVectors(a, b.val);
          else {
            error("Inconsistent dimensions in mul(a,B): size(a) = [" + sa[0] + "," + sa[1] + "], size(B) = [" + sb[0] + "," + sb[1] + "]");
            return undefined;
          }
          break;
        case "spmatrix":
          if (b.m === 1)
            return outerprodVectors(a, fullMatrix(b).val);
          else {
            error("Inconsistent dimensions in mul(a,B): size(a) = [" + sa[0] + "," + sa[1] + "], size(B) = [" + sb[0] + "," + sb[1] + "]");
            return undefined;
          }
          break;
        case "ComplexMatrix":
          if (b.m === 1)
            return transpose(outerprodComplexVectorVector(new ComplexVector(b.re, b.im, true), a, b.val));
          else {
            error("Inconsistent dimensions in mul(a,B): size(a) = [" + sa[0] + "," + sa[1] + "], size(B) = [" + sb[0] + "," + sb[1] + "]");
            return undefined;
          }
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "spvector":
      switch (type(b)) {
        case "number":
          return mulScalarspVector(b, a);
          break;
        case "vector":
          if (a.length !== b.length) {
            error("Error in mul(a,b) (dot product): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return dotspVectorVector(a, b);
          break;
        case "spvector":
          if (a.length !== b.length) {
            error("Error in mul(a,b) (dot product): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return spdot(b, a);
          break;
        case "matrix":
          if (b.m === 1)
            return outerprodspVectorVector(a, b.val);
          else {
            error("Inconsistent dimensions in mul(a,B): size(a) = [" + sa[0] + "," + sa[1] + "], size(B) = [" + sb[0] + "," + sb[1] + "]");
            return undefined;
          }
          break;
        case "spmatrix":
          if (b.m === 1)
            return outerprodspVectorVector(a, fullMatrix(b).val);
          else {
            error("Inconsistent dimensions in mul(a,B): size(a) = [" + sa[0] + "," + sa[1] + "], size(B) = [" + sb[0] + "," + sb[1] + "]");
            return undefined;
          }
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "ComplexVector":
      switch (type(b)) {
        case "number":
          return mulScalarComplexVector(b, a);
          break;
        case "Complex":
          return mulComplexComplexVector(b, a);
          break;
        case "vector":
          if (a.length !== b.length) {
            error("Error in mul(a,b) (dot product): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return dotComplexVectorVector(a, b);
          break;
        case "spvector":
          if (a.length !== b.length) {
            error("Error in mul(a,b) (dot product): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return dotComplexVectorspVector(a, b);
          break;
        case "matrix":
          if (b.m === 1)
            return outerprodComplexVectorVector(a, b.val);
          else {
            error("Inconsistent dimensions in mul(a,B): size(a) = [" + sa[0] + "," + sa[1] + "], size(B) = [" + sb[0] + "," + sb[1] + "]");
            return undefined;
          }
          break;
        case "spmatrix":
          if (b.m === 1)
            return outerprodComplexVectorVector(a, fullMatrix(b).val);
          else {
            error("Inconsistent dimensions in mul(a,B): size(a) = [" + sa[0] + "," + sa[1] + "], size(B) = [" + sb[0] + "," + sb[1] + "]");
            return undefined;
          }
          break;
        case "ComplexMatrix":
          if (b.m === 1)
            return outerprodComplexVectors(a, new ComplexVector(b.re, b.im, true));
          else {
            error("Inconsistent dimensions in mul(a,B): size(a) = [" + sa[0] + "," + sa[1] + "], size(B) = [" + sb[0] + "," + sb[1] + "]");
            return undefined;
          }
          break;
        default:
          return undefined;
          break;
      }
      break;

    case "matrix":
      switch (type(b)) {
        case "number":
          return mulScalarMatrix(b, a);
          break;
        case "Complex":
          return mulComplexMatrix(b, a);
          break;
        case "vector":
          if (a.m === 1) {
            // dot product with explicit transpose
            if (a.val.length !== b.length) {
              error("Error in mul(a',b): a.length = " + a.val.length + " !== " + b.length + " =  b.length.");
              return undefined;
            }
            return dot(a.val, b);
          } else {
            if (a.n !== b.length) {
              error("Error in mul(A,b): A.n = " + a.n + " !== " + b.length + " = b.length.");
              return undefined;
            }
            return mulMatrixVector(a, b);
          }
          break;
        case "spvector":
          if (a.m === 1) {
            // dot product with explicit transpose
            if (a.val.length !== b.length) {
              error("Error in mul(a',b): a.length = " + a.val.length + " !== " + b.length + " =  b.length.");
              return undefined;
            }
            return dotspVectorVector(b, a.val);
          } else {
            if (a.n !== b.length) {
              error("Error in mul(A,b): A.n = " + a.n + " !== " + b.length + " = b.length.");
              return undefined;
            }
            return mulMatrixspVector(a, b);
          }
          break;
        case "ComplexVector":
          if (a.m === 1) {
            // dot product with explicit transpose
            if (a.val.length !== b.length) {
              error("Error in mul(a',b): a.length = " + a.val.length + " !== " + b.length + " =  b.length.");
              return undefined;
            }
            return dotComplexVectorVector(b, a.val);
          } else {
            if (a.n !== b.length) {
              error("Error in mul(A,b): A.n = " + a.n + " !== " + b.length + " = b.length.");
              return undefined;
            }
            return mulMatrixComplexVector(a, b);
          }
          break;
        case "matrix":
          if (a.n !== b.m) {
            error("Error in mul(A,B): A.n = " + a.n + " !== " + b.m + " = B.m.");
            return undefined;
          }
          return mulMatrixMatrix(a, b);
          break;
        case "spmatrix":
          if (a.n !== b.m) {
            error("Error in mul(A,B): A.n = " + a.n + " !== " + b.m + " = B.m.");
            return undefined;
          }
          return mulMatrixspMatrix(a, b);
          break;
        case "ComplexMatrix":
          if (a.n !== b.m) {
            error("Error in mul(A,B): A.n = " + a.n + " !== " + b.m + " = B.m.");
            return undefined;
          }
          return transpose(mulComplexMatrixMatrix(transpose(b), transpose(a)));
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "spmatrix":
      switch (type(b)) {
        case "number":
          return mulScalarspMatrix(b, a);
          break;
        case "vector":
          if (a.m === 1) {
            // dot product with explicit transpose
            if (a.n !== b.length) {
              error("Error in mul(a',b): a.length = " + a.val.length + " !== " + b.length + " =  b.length.");
              return undefined;
            }
            return dot(fullMatrix(a).val, b);
          } else {
            if (a.n !== b.length) {
              error("Error in mul(A,b): A.n = " + a.n + " !== " + b.length + " = b.length.");
              return undefined;
            }
            return mulspMatrixVector(a, b);
          }
          break;
        case "spvector":
          if (a.m === 1) {
            // dot product with explicit transpose
            if (a.n !== b.length) {
              error("Error in mul(a',b): a.length = " + a.val.length + " !== " + b.length + " =  b.length.");
              return undefined;
            }
            return dotspVectorVector(b, fullMatrix(a).val);
          } else {
            if (a.n !== b.length) {
              error("Error in mul(A,b): A.n = " + a.n + " !== " + b.length + " = b.length.");
              return undefined;
            }
            return mulspMatrixspVector(a, b);
          }
          break;
        case "matrix":
          if (a.n !== b.m) {
            error("Error in mul(A,B): A.n = " + a.n + " !== " + b.m + " = B.m.");
            return undefined;
          }
          return mulspMatrixMatrix(a, b);
          break;
        case "spmatrix":
          if (a.n !== b.m) {
            error("Error in mul(A,B): A.n = " + a.n + " !== " + b.m + " = B.m.");
            return undefined;
          }
          return mulspMatrixspMatrix(a, b);
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "ComplexMatrix":
      switch (type(b)) {
        case "number":
          return mulScalarComplexMatrix(b, a);
          break;
        case "Complex":
          return mulComplexComplexMatrix(b, a);
          break;
        case "vector":
          if (a.m === 1) {
            // dot product with explicit transpose
            if (a.val.length !== b.length) {
              error("Error in mul(a',b): a.length = " + a.val.length + " !== " + b.length + " =  b.length.");
              return undefined;
            }
            return dotComplexVectorVector(new ComplexVector(a.re, a.im, true), b);
          } else {
            if (a.n !== b.length) {
              error("Error in mul(A,b): A.n = " + a.n + " !== " + b.length + " = b.length.");
              return undefined;
            }
            return mulComplexMatrixVector(a, b);
          }
          break;
        case "spvector":
          if (a.m === 1) {
            // dot product with explicit transpose
            if (a.val.length !== b.length) {
              error("Error in mul(a',b): a.length = " + a.val.length + " !== " + b.length + " =  b.length.");
              return undefined;
            }
            return dotComplexVectorspVector(new ComplexVector(a.re, a.im, true), b);
          } else {
            if (a.n !== b.length) {
              error("Error in mul(A,b): A.n = " + a.n + " !== " + b.length + " = b.length.");
              return undefined;
            }
            return mulComplexMatrixspVector(a, b);
          }
          break;
        case "ComplexVector":
          if (a.m === 1) {
            // dot product with explicit transpose
            if (a.val.length !== b.length) {
              error("Error in mul(a',b): a.length = " + a.val.length + " !== " + b.length + " =  b.length.");
              return undefined;
            }
            return dotComplexVectors(new ComplexVector(a.re, a.im, true), b);
          } else {
            if (a.n !== b.length) {
              error("Error in mul(A,b): A.n = " + a.n + " !== " + b.length + " = b.length.");
              return undefined;
            }
            return mulComplexMatrixComplexVector(a, b);
          }
          break;
        case "matrix":
          if (a.n !== b.m) {
            error("Error in mul(A,B): A.n = " + a.n + " !== " + b.m + " = B.m.");
            return undefined;
          }
          return mulComplexMatrixMatrix(a, b);
          break;
        case "spmatrix":
          if (a.n !== b.m) {
            error("Error in mul(A,B): A.n = " + a.n + " !== " + b.m + " = B.m.");
            return undefined;
          }
          return mulComplexMatrixspMatrix(a, b);
          break;
        case "ComplexMatrix":
          if (a.n !== b.m) {
            error("Error in mul(A,B): A.n = " + a.n + " !== " + b.m + " = B.m.");
            return undefined;
          }
          return mulComplexMatrices(a, b);
          break;
        default:
          return undefined;
          break;
      }
      break;
    default:
      return undefined;
      break;
  }
}

/**
 * @param {number}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function mulScalarVector(scalar, vec) {
  let i;
  const n = vec.length;
  let res = new Float64Array(vec);
  for (i = 0; i < n; i++) {
    res[i] *= scalar;
  }
  return res;
}

/**
 * @param {number}
 * @param {Matrix}
 * @return {Matrix}
 */
export function mulScalarMatrix(scalar, A) {
  let res = new Matrix(A.m, A.n, mulScalarVector(scalar, A.val), true);

  return res;
}

/**
 * @param {Float64Array}
 * @param {Float64Array}
 * @return {number}
 */
export function dot(a, b) {
  const n = a.length;
  let i;
  let res = 0;
  for (i = 0; i < n; i++) {
    res += a[i]*b[i];
  }
  return res;
}

/**
 * @param {Matrix}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function mulMatrixVector(A, b) {
  const m = A.length;
  let c = new Float64Array(m);
  let r = 0;
  for (let i = 0; i < m; i++) {
    c[i] = dot(A.val.subarray(r, r + A.n), b);
    r += A.n;
  }

  return c;
}

/**
 * @param {Matrix}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function mulMatrixTransVector(A, b) {
  const m = A.length;
  const n = A.n;
  let c = new Float64Array(n);
  let rj = 0;
  for (let j = 0; j < m; j++) {
    let bj = b[j];
    for (let i = 0; i < n; i++) {
      c[i] += A.val[rj + i]*bj;
    }
    rj += A.n;
  }
  return c;
}

/**
 * @param {Matrix}
 * @param {Matrix}
 * @return {Matrix}
 */
export function mulMatrixMatrix(A, B) {
  const m = A.length;
  const n = B.n;
  const n2 = B.length;

  let Av = A.val;
  let Bv = B.val;

  let C = new Float64Array(m*n);
  let aik;
  let Aik = 0;
  let Ci = 0;
  for (let i = 0; i < m; i++) {
    let bj = 0;
    for (let k = 0; k < n2; k++) {
      aik = Av[Aik];
      for (let j = 0; j < n; j++) {
        C[Ci + j] += aik*Bv[bj];
        bj++;
      }
      Aik++;
    }
    Ci += n;
  }
  return new Matrix(m, n, C, true);
}

/**
 * @param {Float64Array}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function entrywisemulVector(a, b) {
  let i;
  const n = a.length;
  let res = new Float64Array(n);
  for (i = 0; i < n; i++) {
    res[i] = a[i]*b[i];
  }
  return res;
}

/**
 * @param {Matrix}
 * @param {Matrix}
 * @return {Matrix}
 */
export function entrywisemulMatrix(A, B) {
  let res = new Matrix(A.m, A.n, entrywisemulVector(A.val, B.val), true);
  return res;
}


export function entrywisemul(a, b) {
  let sa = size(a);
  let sb = size(b);
  if (typeof (a) !== "number" && sa[0] === 1 && sa[1] === 1)
    a = get(a, 0, 0);
  if (typeof (b) !== "number" && sb[0] === 1 && sb[1] === 1)
    b = get(b, 0, 0);

  switch (type(a)) {
    case "number":
      switch (type(b)) {
        case "number":
          return a*b;
          break;
        case "Complex":
          return mulComplexReal(b, a);
          break;
        case "vector":
          return mulScalarVector(a, b);
          break;
        case "spvector":
          return mulScalarspVector(a, b);
          break;
        case "ComplexVector":
          return mulScalarComplexVector(b, a);
          break;
        case "matrix":
          return mulScalarMatrix(a, b);
          break;
        case "spmatrix":
          return mulScalarspMatrix(a, b);
          break;
        case "ComplexMatrix":
          return mulScalarComplexMatrix(b, a);
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "vector":
      switch (type(b)) {
        case "number":
          return mulScalarVector(b, a);
          break;
        case "Complex":
          return mulComplexVector(b, a);
          break;
        case "vector":
          if (a.length !== b.length) {
            error("Error in entrywisemul(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return entrywisemulVector(a, b);
          break;
        case "ComplexVector":
          if (a.length !== b.length) {
            error("Error in entrywisemul(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return entrywisemulComplexVectorVector(b, a);
          break;
        case "spvector":
          if (a.length !== b.length) {
            error("Error in entrywisemul(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return entrywisemulspVectorVector(b, a);
          break;
        case "matrix":
        case "spmatrix":
        case "ComplexMatrix":
          error("Error in entrywisemul(a,B): a is a vector and B is a matrix.");
          return undefined;
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "spvector":
      switch (type(b)) {
        case "number":
          return mulScalarspVector(b, a);
          break;
        case "vector":
          if (a.length !== b.length) {
            error("Error in entrywisemul(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return entrywisemulspVectorVector(a, b);
          break;
        case "spvector":
          if (a.length !== b.length) {
            error("Error in entrywisemul(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return entrywisemulspVectors(a, b);
          break;
        case "matrix":
          error("Error in entrywisemul(a,B): a is a vector and B is a Matrix.");
          return undefined;
          break;
        case "spmatrix":
          error("Error in entrywisemul(a,B): a is a vector and B is a Matrix.");
          return undefined;
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "matrix":
      switch (type(b)) {
        case "number":
          return mulScalarMatrix(b, a);
          break;
        case "Complex":
          return mulComplexMatrix(b, a);
          break;
        case "vector":
        case "spvector":
        case "ComplexVector":
          error("Error in entrywisemul(A,b): A is a Matrix and b is a vector.");
          return undefined;
          break;
        case "matrix":
          if (a.m !== b.m || a.n !== b.n) {
            error("Error in entrywisemul(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
            return undefined;
          }
          return entrywisemulMatrix(a, b);
          break;
        case "spmatrix":
          if (a.m !== b.m || a.n !== b.n) {
            error("Error in entrywisemul(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
            return undefined;
          }
          return entrywisemulspMatrixMatrix(b, a);
          break;
        case "ComplexMatrix":
          if (a.m !== b.m || a.n !== b.n) {
            error("Error in entrywisemul(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
            return undefined;
          }
          return entrywisemulComplexMatrixMatrix(b, a);
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "spmatrix":
      switch (type(b)) {
        case "number":
          return mulScalarspMatrix(b, a);
          break;
        case "vector":
          error("Error in entrywisemul(A,b): A is a Matrix and b is a vector.");
          return undefined;
          break;
        case "spvector":
          error("Error in entrywisemul(A,b): A is a Matrix and b is a vector.");
          return undefined;
          break;
        case "matrix":
          if (a.m !== b.m || a.n !== b.n) {
            error("Error in entrywisemul(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
            return undefined;
          }
          return entrywisemulspMatrixMatrix(a, b);
          break;
        case "spmatrix":
          if (a.m !== b.m || a.n !== b.n) {
            error("Error in entrywisemul(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
            return undefined;
          }
          return entrywisemulspMatrices(a, b);
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "ComplexVector":
      switch (type(b)) {
        case "number":
          return mulScalarComplexVector(b, a);
          break;
        case "Complex":
          return mulComplexComplexVector(b, a);
          break;
        case "vector":
          if (a.length !== b.length) {
            error("Error in entrywisemul(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return entrywisemulComplexVectorVector(a, b);
          break;
        case "ComplexVector":
          if (a.length !== b.length) {
            error("Error in entrywisemul(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return entrywisemulComplexVectors(a, b);
          break;
        case "spvector":
          if (a.length !== b.length) {
            error("Error in entrywisemul(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return entrywisemulComplexVectorspVector(a, b);
          break;
        case "matrix":
        case "spmatrix":
        case "ComplexMatrix":
          error("Error in entrywisemul(a,B): a is a vector and B is a matrix.");
          return undefined;
          break;
        default:
          return undefined;
          break;
      }
      break;
    case "ComplexMatrix":
      switch (type(b)) {
        case "number":
          return mulScalarComplexMatrix(b, a);
          break;
        case "Complex":
          return mulComplexComplexMatrix(b, a);
          break;
        case "vector":
        case "spvector":
        case "ComplexVector":
          error("Error in entrywisemul(A,b): A is a Matrix and b is a vector.");
          return undefined;
          break;
        case "matrix":
          if (a.m !== b.m || a.n !== b.n) {
            error("Error in entrywisemul(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
            return undefined;
          }
          return entrywisemulComplexMatrixMatrix(a, b);
          break;
        case "spmatrix":
          if (a.m !== b.m || a.n !== b.n) {
            error("Error in entrywisemul(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
            return undefined;
          }
          return entrywisemulComplexMatrixspMatrix(a, b);
          break;
        case "ComplexMatrix":
          if (a.m !== b.m || a.n !== b.n) {
            error("Error in entrywisemul(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
            return undefined;
          }
          return entrywisemulComplexMatrices(a, b);
          break;
        default:
          return undefined;
          break;
      }
      break;
    default:
      return undefined;
      break;
  }
}


/** SAXPY : y = y + ax
 * @param {number}
 * @param {Float64Array}
 * @param {Float64Array}
 */
export function saxpy(a, x, y) {
  const n = y.length;
  for (let i = 0; i < n; i++) {
    y[i] += a*x[i];
  }
}

/** GAXPY : y = y + Ax
 * @param {Matrix}
 * @param {Float64Array}
 * @param {Float64Array}
 */
export function gaxpy(A, x, y) {
  const m = A.m;
  const n = A.n;
  let r = 0;
  for (let i = 0; i < m; i++) {
    y[i] += dot(A.val.subarray(r, r + n), x);
    r += n;
  }
}

/**
 * @param {Float64Array}
 * @param {number}
 * @return {Float64Array}
 */
export function divVectorScalar(a, b) {
  let i;
  const n = a.length;
  let res = new Float64Array(a);
  for (i = 0; i < n; i++) {
    res[i] /= b;
  }
  return res;
}

/**
 * @param {number}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function divScalarVector(a, b) {
  let i;
  const n = b.length;
  let res = new Float64Array(n);
  for (i = 0; i < n; i++) {
    res[i] = a/b[i];
  }
  return res;
}

/**
 * @param {Float64Array}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function divVectors(a, b) {
  let i;
  const n = a.length;
  let res = new Float64Array(a);
  for (i = 0; i < n; i++) {
    res[i] /= b[i];
  }
  return res;
}

/**
 * @param {Matrix}
 * @param {number}
 * @return {Matrix}
 */
export function divMatrixScalar(A, b) {
  let res = new Matrix(A.m, A.n, divVectorScalar(A.val, b), true);
  return res;
}

/**
 * @param {number}
 * @param {Matrix}
 * @return {Matrix}
 */
export function divScalarMatrix(a, B) {
  let res = new Matrix(B.m, B.n, divScalarVector(a, B.val), true);
  return res;
}

/**
 * @param {Matrix}
 * @param {Matrix}
 * @return {Matrix}
 */
export function divMatrices(A, B) {
  let res = new Matrix(A.m, A.n, divVectors(A.val, B.val), true);
  return res;
}

export function entrywisediv(a, b) {
  let ta = type(a);
  let tb = type(b);

  switch (ta) {
    case "number":
      switch (tb) {
        case "number":
          return a/b;
          break;
        case "vector":
          return divScalarVector(a, b);
          break;
        case "matrix":
          return divScalarMatrix(a, b);
          break;
        case "spvector":
          return divScalarspVector(a, b);
          break;
        case "spmatrix":
          return divScalarspMatrix(a, b);
          break;
        default:
          error("Error in entrywisediv(a,b): b must be a number, a vector or a matrix.");
          return undefined;
      }
      break;
    case "vector":
      switch (tb) {
        case "number":
          return divVectorScalar(a, b);
          break;
        case "vector":
          if (a.length !== b.length) {
            error("Error in entrywisediv(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return divVectors(a, b);
          break;
        case "spvector":
          error("Error in entrywisediv(a,b): b is a sparse vector with zeros.");
          break;
        default:
          error("Error in entrywisediv(a,B): a is a vector and B is a " + tb + ".");
          return undefined;
      }
      break;
    case "spvector":
      switch (tb) {
        case "number":
          return mulScalarspVector(1/b, a);
          break;
        case "vector":
          if (a.length !== b.length) {
            error("Error in entrywisediv(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
            return undefined;
          }
          return divVectorspVector(a, b);
          break;
        case "spvector":
          error("Error in entrywisediv(a,b): b is a sparse vector with zeros.");
          return undefined;
          break;
        default:
          error("Error in entrywisediv(a,B): a is a vector and B is a " + tb + ".");
          return undefined;
      }
      break;
    case "matrix":
      switch (tb) {
        case "number":
          return divMatrixScalar(a, b);
          break;
        case "matrix":
          if (a.m !== b.m || a.n !== b.n) {
            error("Error in entrywisediv(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
            return undefined;
          }
          return divMatrices(a, b);
          break;
        case "spmatrix":
          error("Error in entrywisediv(A,B): B is a sparse matrix with zeros.");
          return undefined;
          break;
        default:
          error("Error in entrywisediv(A,b): a is a matrix and B is a " + tb + ".");
          return undefined;
      }
    case "spmatrix":
      switch (tb) {
        case "number":
          return mulScalarspMatrix(1/b, a);
          break;
        case "matrix":
          if (a.m !== b.m || a.n !== b.n) {
            error("Error in entrywisediv(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
            return undefined;
          }
          return divMatrixspMatrix(a, b);
          break;
        case "spmatrix":
          error("Error in entrywisediv(A,B): B is a sparse matrix with zeros.");
          return undefined;
          break;
        default:
          error("Error in entrywisediv(A,b): a is a matrix and B is a " + tb + ".");
          return undefined;
      }
      break;
    default:
      error("Error in entrywisediv(a,b): a must be a number, a vector or a matrix.");
      return undefined;
      break;
  }
}

export function outerprodVectors(a, b, scalar) {
  let i;
  let j;
  let ui;
  const m = a.length;
  const n = b.length;
  let res = new Matrix(m, n);
  if (arguments.length === 3) {
    for (i = 0; i < m; i++) {
      res.val.set(mulScalarVector(scalar*a[i], b), i*n);
    }
  } else {
    for (i = 0; i < m; i++) {
      res.val.set(mulScalarVector(a[i], b), i*n);
    }
  }
  return res;
}

export function outerprod(u, v, scalar) {
  // outer product of two vectors : res = scalar * u * v^T

  if (typeof (u) === "number") {
    if (typeof (v) === "number") {
      if (arguments.length === 2)
        return u*v;
      else
        return u*v*scalar;
    } else {
      if (arguments.length === 2)
        return new Matrix(1, v.length, mulScalarVector(u, v), true);
      else
        return new Matrix(1, v.length, mulScalarVector(u*scalar, v), true);
    }
  }
  if (u.length === 1) {
    if (typeof (v) === "number") {
      if (arguments.length === 2)
        return u[0]*v;
      else
        return u[0]*v*scalar;
    } else {
      if (arguments.length === 2)
        return new Matrix(1, v.length, mulScalarVector(u[0], v), true);
      else
        return new Matrix(1, v.length, mulScalarVector(u[0]*scalar, v), true);
    }
  }
  if (typeof (v) === "number") {
    if (arguments.length === 2)
      return mulScalarVector(v, u);
    else
      return mulScalarVector(scalar*v, u);
  }
  if (v.length === 1) {
    if (arguments.length === 2)
      return mulScalarVector(v[0], u);
    else
      return mulScalarVector(scalar*v[0], u);
  }

  if (arguments.length === 2)
    return outerprodVectors(u, v);
  else
    return outerprodVectors(u, v, scalar);
}

/**
 * @param {number}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function addScalarVector(scalar, vec) {
  const n = vec.length;
  let res = new Float64Array(vec);
  for (let i = 0; i < n; i++) {
    res[i] += scalar;
  }

  return res;
}

/**
 * @param {number}
 * @param {Matrix}
 * @return {Matrix}
 */
export function addScalarMatrix(a, B) {
  return new Matrix(B.m, B.n, addScalarVector(a, B.val), true);
}

/**
 * @param {Float64Array}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function addVectors(a, b) {
  const n = a.length;
  let c = new Float64Array(a);
  for (let i = 0; i < n; i++) {
    c[i] += b[i];
  }
  return c;
}

/**
 * @param {Matrix}
 * @param {Matrix}
 * @return {Matrix}
 */
export function addMatrices(A, B) {
  return new Matrix(A.m, A.n, addVectors(A.val, B.val), true);
}

export function add(a, b) {

  const ta = type(a);
  const tb = type(b);
  if (ta === "number" && tb === "number" || ta === "string" || tb === "string")
    return a + b;
  else if (ta === "number") {
    switch (tb) {
      case "Complex":
        return addComplexReal(b, a);
        break;
      case "vector":
        return addScalarVector(a, b);
        break;
      case "matrix":
        return addScalarMatrix(a, b);
        break;
      case "spvector":
        return addScalarspVector(a, b);
        break;
      case "spmatrix":
        return addScalarspMatrix(a, b);
        break;
      case "ComplexVector":
        return addScalarComplexVector(a, b);
        break;
      case "ComplexMatrix":
        return addScalarComplexMatrix(a, b);
        break;
      default:
        return undefined;
        break;
    }
  } else if (tb === "number") {
    switch (ta) {
      case "Complex":
        return addComplexReal(a, b);
        break;
      case "vector":
        return addScalarVector(b, a);
        break;
      case "matrix":
        return addScalarMatrix(b, a);
        break;
      case "spvector":
        return addScalarspVector(b, a);
        break;
      case "spmatrix":
        return addScalarspMatrix(b, a);
        break;
      case "ComplexVector":
        return addScalarComplexVector(b, a);
        break;
      case "ComplexMatrix":
        return addScalarComplexMatrix(b, a);
        break;
      default:
        return undefined;
        break;
    }
  } else if (ta === "vector") {
    switch (tb) {
      case "vector":
        // vector addition
        if (a.length !== b.length) {
          error("Error in add(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return addVectors(a, b);
        break;
      case "spvector":
        if (a.length !== b.length) {
          error("Error in add(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return addVectorspVector(a, b);
        break;
      case "ComplexVector":
        if (a.length !== b.length) {
          error("Error in add(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return addComplexVectorVector(b, a);
        break;
      case "matrix":
      case "spmatrix":
      default:
        error("Error in add(a,B): a is a vector and B is a " + tb + ".");
        return undefined;
        break;
    }
  } else if (ta === "matrix") {
    switch (tb) {
      case "matrix":
        // Matrix addition
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in add(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return addMatrices(a, b);
        break;
      case "spmatrix":
        // Matrix addition
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in add(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return addMatrixspMatrix(a, b);
        break;
      case "ComplexMatrix":
        // Matrix addition
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in add(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return addComplexMatrixMatrix(b, a);
        break;
      case "vector":
      case "spvector":
      default:
        error("Error in add(A,b): a is a matrix and B is a " + tb + ".");
        return undefined;
        break;
    }
  } else if (ta === "spvector") {
    switch (tb) {
      case "vector":
        // vector addition
        if (a.length !== b.length) {
          error("Error in add(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return addVectorspVector(b, a);
        break;
      case "spvector":
        if (a.length !== b.length) {
          error("Error in add(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return addspVectors(a, b);
        break;
      case "matrix":
      case "spmatrix":
      default:
        error("Error in add(a,B): a is a sparse vector and B is a " + tb + ".");
        return undefined;
        break;
    }
  } else if (ta === "spmatrix") {
    switch (tb) {
      case "matrix":
        // Matrix addition
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in add(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return addMatrixspMatrix(b, a);
        break;
      case "spmatrix":
        // Matrix addition
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in add(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return addspMatrices(a, b);
        break;
      case "vector":
      case "spvector":
      default:
        error("Error in add(A,b): a is a sparse matrix and B is a " + tb + ".");
        return undefined;
        break;
    }
  } else if (ta === "ComplexVector") {
    switch (tb) {
      case "vector":
        // vector addition
        if (a.length !== b.length) {
          error("Error in add(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return addComplexVectorVector(a, b);
        break;
      case "spvector":
        if (a.length !== b.length) {
          error("Error in add(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return addComplexVectorspVector(a, b);
        break;
      case "ComplexVector":
        if (a.length !== b.length) {
          error("Error in add(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return addComplexVectors(b, a);
        break;
      case "matrix":
      case "spmatrix":
      default:
        error("Error in add(a,B): a is a vector and B is a " + tb + ".");
        return undefined;
        break;
    }
  } else if (ta === "ComplexMatrix") {
    switch (tb) {
      case "matrix":
        // Matrix addition
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in add(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return addComplexMatrixMatrix(a, b);
        break;
      case "spmatrix":
        // Matrix addition
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in add(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return addComplexMatrixspMatrix(a, b);
        break;
      case "ComplexMatrix":
        // Matrix addition
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in add(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return addComplexMatrices(a, b);
        break;
      case "vector":
      case "spvector":
      default:
        error("Error in add(A,b): a is a matrix and B is a " + tb + ".");
        return undefined;
        break;
    }
  } else
    return undefined;
}

/**
 * @param {number}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function subScalarVector(scalar, vec) {
  const n = vec.length;
  let res = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    res[i] = scalar - vec[i];
  }

  return res;
}

/**
 * @param {Float64Array}
 * @param {number}
 * @return {Float64Array}
 */
export function subVectorScalar(vec, scalar) {
  const n = vec.length;
  let res = new Float64Array(vec);
  for (let i = 0; i < n; i++) {
    res[i] -= scalar;
  }

  return res;
}

/**
 * @param {number}
 * @param {Matrix}
 * @return {Matrix}
 */
export function subScalarMatrix(a, B) {
  return new Matrix(B.m, B.n, subScalarVector(a, B.val), true);
}

/**
 * @param {Matrix}
 * @param {number}
 * @return {Matrix}
 */
export function subMatrixScalar(B, a) {
  return new Matrix(B.m, B.n, subVectorScalar(B.val, a), true);
}

/**
 * @param {Float64Array}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function subVectors(a, b) {
  const n = a.length;
  let c = new Float64Array(a);
  for (let i = 0; i < n; i++) {
    c[i] -= b[i];
  }
  return c;
}

/**
 * @param {Matrix}
 * @param {Matrix}
 * @return {Matrix}
 */
export function subMatrices(A, B) {
  return new Matrix(A.m, A.n, subVectors(A.val, B.val), true);
}

export function sub(a, b) {

  const ta = type(a);
  const tb = type(b);
  if (ta === "number" && tb === "number")
    return a - b;
  else if (ta === "number") {
    switch (tb) {
      case "Complex":
        return addComplexReal(minusComplex(b), a);
        break;
      case "vector":
        return subScalarVector(a, b);
        break;
      case "matrix":
        return subScalarMatrix(a, b);
        break;
      case "spvector":
        return subScalarspVector(a, b);
        break;
      case "spmatrix":
        return subScalarspMatrix(a, b);
        break;
      default:
        return undefined;
        break;
    }
  } else if (tb === "number") {
    switch (ta) {
      case "Complex":
        return addComplexReal(b, -a);
        break;
      case "vector":
        return subVectorScalar(a, b);
        break;
      case "matrix":
        return subMatrixScalar(a, b);
        break;
      case "spvector":
        return addScalarspVector(-b, a);
        break;
      case "spmatrix":
        return addScalarspMatrix(-b, a);
        break;
      default:
        return undefined;
        break;
    }
  } else if (ta === "vector") {
    switch (tb) {
      case "vector":
        // vector substraction
        if (a.length !== b.length) {
          error("Error in sub(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return subVectors(a, b);
        break;
      case "spvector":
        // vector substraction
        if (a.length !== b.length) {
          error("Error in sub(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return subVectorspVector(a, b);
        break;
      case "matrix":
      case "spmatrix":
      default:
        error("Error in sub(a,B): a is a vector and B is a " + tb + ".");
        return undefined;
        break;
    }
  } else if (ta === "matrix") {
    switch (tb) {
      case "matrix":
        // Matrix sub
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in sub(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return subMatrices(a, b);
        break;
      case "spmatrix":
        // Matrix addition
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in sub(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return subMatrixspMatrix(a, b);
        break;
      case "vector":
      case "spvector":
      default:
        error("Error in sub(A,b): A is a matrix and b is a " + tb + ".");
        return undefined;
        break;
    }
  } else if (ta === "spvector") {
    switch (tb) {
      case "vector":
        if (a.length !== b.length) {
          error("Error in sub(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return subspVectorVector(a, b);
        break;
      case "spvector":
        if (a.length !== b.length) {
          error("Error in sub(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
          return undefined;
        }
        return subspVectors(a, b);
        break;
      case "matrix":
      case "spmatrix":
      default:
        error("Error in sub(a,B): a is a sparse vector and B is a " + tb + ".");
        return undefined;
        break;
    }
  } else if (ta === "spmatrix") {
    switch (tb) {
      case "matrix":
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in sub(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return subspMatrixMatrix(a, b);
        break;
      case "spmatrix":
        if (a.m !== b.m || a.n !== b.n) {
          error("Error in sub(A,B): size(A) = [" + a.m + "," + a.n + "] !== [" + b.m + "," + b.n + "] = size(B).");
          return undefined;
        }
        return subspMatrices(a, b);
        break;
      case "vector":
      case "spvector":
      default:
        error("Error in sub(A,b): a is a sparse matrix and B is a " + tb + ".");
        return undefined;
        break;
    }
  } else
    return undefined;
}

export function pow(a, b) {
  let i;
  const ta = type(a);
  const tb = type(b);

  if (ta === "number" && tb === "number")
    return Math.pow(a, b);
  else if (ta === "number") {
    if (tb === "vector") {
      let c = zeros(b.length);
      if (!isZero(a)) {
        for (i = 0; i < b.length; i++) {
          c[i] = Math.pow(a, b[i]);
        }
      }
      return c;
    } else {
      let c = new Matrix(b.m, b.n, pow(a, b.val), true);
      return c;
    }
  } else if (tb === "number") {
    if (ta === "vector") {
      let c = zeros(a.length);
      for (i = 0; i < a.length; i++) {
        c[i] = Math.pow(a[i], b);
      }
      return c;
    } else {
      let c = new Matrix(a.m, a.n, pow(a.val, b), true);
      return c;
    }
  } else if (ta === "vector") {
    if (tb === "vector") {
      // entry-wise power
      if (a.length !== b.length) {
        error("Error in pow(a,b): a.length = " + a.length + " !== " + b.length + " = b.length.");
        return undefined;
      }
      let c = zeros(a.length);
      for (i = 0; i < a.length; i++) {
        c[i] = Math.pow(a[i], b[i]);
      }
      return c;
    } else {
      // vector + matrix
      return "undefined";
    }
  } else {
    if (tb === "vector") {
      // matrix + vector
      return "undefined";
    } else {
      // entry-wise power
      let c = new Matrix(a.m, a.n, pow(a.val, b.val), true);
      return c;
    }
  }
}

export function minus(x) {

  switch (type(x)) {
    case "number":
      return -x;
    case "vector":
      return minusVector(x);
    case "spvector":
      return new spVector(x.length, minusVector(x.val), x.ind);
    case "ComplexVector":
      return minusComplexVector(x);
    case "matrix":
      return new Matrix(x.m, x.n, minusVector(x.val), true);
    case "spmatrix":
      return new spMatrix(x.m, x.n, minusVector(x.val), x.cols, x.rows);
    case "ComplexMatrix":
      return minusComplexMatrix(x);
    default:
      return undefined;
  }
}

/**
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function minusVector(x) {
  let res = new Float64Array(x.length);
  for (let i = 0; i < x.length; i++) {
    res[i] = -x[i];
  }
  return res;
}

/**
 * @param {Matrix}
 * @return {Matrix}
 */
export function minusMatrix(x) {
  return new Matrix(x.m, x.n, minusVector(x.val), true);
}

// minimum

/**
 * @param {Float64Array}
 * @return {number}
 */
export function minVector(a) {
  const n = a.length;
  let res = a[0];
  for (let i = 1; i < n; i++) {
    if (a[i] < res)
      res = a[i];
  }
  return res;
}

/**
 * @param {Matrix}
 * @return {number}
 */
export function minMatrix(A) {
  return minVector(A.val);
}

/**
 * @param {Float64Array}
 * @param {number}
 * @return {Float64Array}
 */
export function minVectorScalar(vec, scalar) {
  let n = vec.length;
  let res = new Float64Array(vec);
  for (let i = 0; i < n; i++) {
    if (scalar < vec[i])
      res[i] = scalar;
  }
  return res;
}

/**
 * @param {Matrix}
 * @param {number}
 * @return {Matrix}
 */
export function minMatrixScalar(A, scalar) {
  return new Matrix(A.m, A.n, minVectorScalar(A.val, scalar), true);
}

/**
 * @param {Matrix}
 * @return {Matrix}
 */
export function minMatrixRows(A) {
  const m = A.m;
  const n = A.n;
  let res = new Float64Array(A.val.subarray(0, n));
  let j;
  let r = n;
  for (let i = 1; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (A.val[r + j] < res[j])
        res[j] = A.val[r + j];
    }
    r += n;
  }
  return new Matrix(1, n, res, true);
}

/**
 * @param {Matrix}
 * @return {Float64Array}
 */
export function minMatrixCols(A) {
  let m = A.m;
  let res = new Float64Array(m);
  let r = 0;
  for (let i = 0; i < m; i++) {
    res[i] = minVector(A.val.subarray(r, r + A.n));
    r += A.n;
  }
  return res;
}

/**
 * @param {Float64Array}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function minVectorVector(a, b) {
  const n = a.length;
  let res = new Float64Array(a);
  for (let i = 0; i < n; i++) {
    if (b[i] < a[i])
      res[i] = b[i];
  }
  return res;
}

/**
 * @param {Matrix}
 * @param {Matrix}
 * @return {Matrix}
 */
export function minMatrixMatrix(A, B) {
  return new Matrix(A.m, A.n, minVectorVector(A.val, B.val), true);
}

export function min(a, b) {
  let ta = type(a);

  if (arguments.length === 1) {
    switch (ta) {
      case "vector": {
        return minVector(a);
      }
      case "spvector": {
        let m = minVector(a.val);
        if (m > 0 && a.val.length < a.length)
          return 0;
        else
          return m;
      }
      case "matrix": {
        return minMatrix(a);
      }
      case "spmatrix": {
        let m = minVector(a.val);

        if (m > 0 && a.val.length < a.m*a.n)
          return 0;
        else
          return m;
      }
      default:
        return a;
    }
  }

  let tb = type(b);
  if (ta === "spvector") {
    a = fullVector(a);
    ta = "vector";
  }
  if (ta === "spmatrix") {
    a = fullMatrix(a);
    ta = "matrix";
  }
  if (tb === "spvector") {
    b = fullVector(b);
    tb = "vector";
  }
  if (tb === "spmatrix") {
    b = fullMatrix(b);
    tb = "matrix";
  }

  if (ta === "number" && tb === "number")
    return Math.min(a, b);
  else if (ta === "number") {
    if (tb === "vector")
      return minVectorScalar(b, a);
    else
      return minMatrixScalar(b, a);
  } else if (tb === "number") {
    if (ta === "vector")
      return minVectorScalar(a, b);
    else {
      // MAtrix , scalar
      if (b === 1)
        return minMatrixRows(a); // return row vector of min of columns
      else if (b === 2)
        return minMatrixCols(a); // return column vector of min of rows
      else
        return minMatrixScalar(a, b);
    }
  } else if (ta === "vector") {
    if (tb === "vector")
      return minVectorVector(a, b);
    else
      return "undefined";
  } else {
    if (tb === "matrix")
      return minMatrixMatrix(a, b);
    else
      return "undefined";
  }
}

// maximum
/**
 * @param {Float64Array}
 * @return {number}
 */
export function maxVector(a) {
  const n = a.length;
  let res = a[0];
  for (let i = 1; i < n; i++) {
    if (a[i] > res)
      res = a[i];
  }
  return res;
}

/**
 * @param {Matrix}
 * @return {number}
 */
export function maxMatrix(A) {
  return maxVector(A.val);
}

/**
 * @param {Float64Array}
 * @param {number}
 * @return {Float64Array}
 */
export function maxVectorScalar(vec, scalar) {
  const n = vec.length;
  let res = new Float64Array(vec);
  for (let i = 0; i < n; i++) {
    if (scalar > vec[i])
      res[i] = scalar;
  }
  return res;
}

/**
 * @param {Matrix}
 * @param {number}
 * @return {Matrix}
 */
export function maxMatrixScalar(A, scalar) {
  return maxVectorScalar(A.val, scalar);
}

/**
 * @param {Matrix}
 * @return {Matrix}
 */
export function maxMatrixRows(A) {
  const m = A.m;
  const n = A.n;
  let res = new Float64Array(A.val.subarray(0, n));
  let j;
  let r = n;
  for (let i = 1; i < m; i++) {
    for (j = 0; j < n; j++) {
      if (A.val[r + j] > res[j])
        res[j] = A.val[r + j];
    }
    r += n;
  }
  return new Matrix(1, n, res, true);
}

/**
 * @param {Matrix}
 * @return {Float64Array}
 */
export function maxMatrixCols(A) {
  const m = A.m;
  let res = new Float64Array(m);
  let r = 0;
  for (let i = 0; i < m; i++) {
    res[i] = maxVector(A.val.subarray(r, r + A.n));
    r += A.n;
  }
  return res;
}

/**
 * @param {Float64Array}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function maxVectorVector(a, b) {
  let n = a.length;
  let res = new Float64Array(a);
  for (let i = 0; i < n; i++) {
    if (b[i] > a[i])
      res[i] = b[i];
  }
  return res;
}

/**
 * @param {Matrix}
 * @param {Matrix}
 * @return {Matrix}
 */
export function maxMatrixMatrix(A, B) {
  return new Matrix(A.m, A.n, maxVectorVector(A.val, B.val), true);
}

export function max(a, b) {
  let ta = type(a);

  if (arguments.length === 1) {
    switch (ta) {
      case "vector":
        return maxVector(a);
      case "spvector": {
        let m = maxVector(a.val);
        if (m < 0 && a.val.length < a.length)
          return 0;
        else
          return m;
      }
      case "matrix":
        return maxMatrix(a);
      case "spmatrix": {
        let m = maxVector(a.val);
        if (m < 0 && a.val.length < a.m*a.n)
          return 0;
        else
          return m;
      }
      default:
        return a;
    }
  }

  let tb = type(b);
  if (ta === "spvector") {
    a = fullVector(a);
    ta = "vector";
  }
  if (ta === "spmatrix") {
    a = fullMatrix(a);
    ta = "matrix";
  }
  if (tb === "spvector") {
    b = fullVector(b);
    tb = "vector";
  }
  if (tb === "spmatrix") {
    b = fullMatrix(b);
    tb = "matrix";
  }

  if (ta === "number" && tb === "number")
    return Math.max(a, b);
  else if (ta === "number") {
    if (tb === "vector")
      return maxVectorScalar(b, a);
    else
      return maxMatrixScalar(b, a);
  } else if (tb === "number") {
    if (ta === "vector")
      return maxVectorScalar(a, b);
    else {
      // MAtrix , scalar
      if (b === 1)
        return maxMatrixRows(a); // return row vector of max of columns
      else if (b === 2)
        return maxMatrixCols(a); // return column vector of max of rows
      else
        return maxMatrixScalar(a, b);
    }
  } else if (ta === "vector") {
    if (tb === "vector")
      return maxVectorVector(a, b);
    else
      return "undefined";
  } else {
    if (tb === "matrix")
      return maxMatrixMatrix(a, b);
    else
      return "undefined";
  }
}

/**
 * @param {Matrix}
 */
export function transposeMatrix(A) {
  let i;
  let j;
  const m = A.m;
  const n = A.n;
  if (m > 1) {
    let res = zeros(n, m);
    let Aj = 0;
    for (j = 0; j < m; j++) {
      let ri = 0;
      for (i = 0; i < n; i++) {
        res.val[ri + j] = A.val[Aj + i];
        ri += m;
      }
      Aj += n;
    }
    return res;
  } else {
    return A.val;
  }
}

/**
 * @param {Float64Array}
 * @return {Matrix}
 */
export function transposeVector(a) {
  return new Matrix(1, a.length, a);
}

export function transpose(A) {
  let i;
  let j;

  switch (type(A)) {
    case "number":
      return A;
      break;
    case "vector": {
      let res = new Matrix(1, A.length, A);
      return res;	// matrix with a single row
    }
    case "spvector":
      return transposespVector(A);
    case "ComplexVector": {
      let res = new ComplexMatrix(1, A.length, conj(A));
      return res;	// matrix with a single row
    }
    case "matrix":
      return transposeMatrix(A);
    case "spmatrix":
      return transposespMatrix(A);
    case "ComplexMatrix":
      return transposeComplexMatrix(A);
    default:
      return undefined;
  }
}

/**
 * @param {Matrix}
 * @return {number}
 */
export function det(A) {
  const n = A.n;
  if (A.m !== n || typeof (A.m) == "undefined")
    return undefined;

  if (n === 2) {
    return A.val[0]*A.val[3] - A.val[1]*A.val[2];
  } else {
    let detA = 0;
    let i, j;
    for (i = 0; i < n; i++) {
      let proddiag = 1;
      for (j = 0; j < n; j++) {
        proddiag *= A.val[((i + j)%n)*n + j];
      }

      detA += proddiag;
    }
    for (i = 0; i < n; i++) {
      let proddiag = 1;
      for (j = 0; j < n; j++) {
        proddiag *= A.val[((i + n - 1 - j)%n)*n + j];
      }

      detA -= proddiag;
    }
  }
  return detA;
}

export function trace(A) {
  if (type(A) === "matrix") {
    let n = A.length;
    if (A.m !== n)
      return "undefined";
    let res = 0;
    for (let i = 0; i < n; i++) {
      res += A.val[i*n + i];
    }
    return res;
  } else {
    return undefined;
  }
}

/**
 * @param {Matrix}
 * @return {Matrix}
 */
export function triu(A) {
  // return the upper triangular part of A
  let i;
  let j;
  const n = A.n;
  const m = A.m;
  let res = zeros(m, n);
  let im = m;
  if (n < m)
    im = n;
  let r = 0;
  for (i = 0; i < im; i++) {
    for (j = i; j < n; j++) {
      res.val[r + j] = A.val[r + j];
    }
    r += n;
  }
  return res;
}

/**
 * @param {Matrix}
 * @return {Matrix}
 */
export function tril(A) {
  // return the lower triangular part of A
  let i;
  let j;
  const n = A.n;
  const m = A.m;
  let res = zeros(m, n);
  let im = m;
  if (n < m)
    im = n;
  let r = 0;
  for (i = 0; i < im; i++) {
    for (j = 0; j <= i; j++) {
      res.val[r + j] = A.val[r + j];
    }
    r += n;
  }
  if (m > im) {
    for (i = im; i < m; i++) {
      for (j = 0; j < n; j++) {
        res.val[r + j] = A.val[r + j];
      }
      r += n;
    }
  }
  return res;
}

/**
 * @param {Matrix}
 * @return {boolean}
 */
export function issymmetric(A) {
  const m = A.m;
  const n = A.n;
  if (m !== n)
    return false;

  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      if (A.val[i*n + j] !== A.val[j*n + i])
        return false;
    }
  }

  return true;
}

/** Concatenate matrices/vectors
 * @param {Array}
 * @param {boolean}
 * @return {Matrix}
 */
export function mat(elems, rowwise) {
  let k;
  let concatWithNumbers = false;
  let elemtypes = new Array(elems.length);
  for (k = 0; k < elems.length; k++) {
    elemtypes[k] = type(elems[k]);
    if (elemtypes[k] === "number")
      concatWithNumbers = true;
  }


  if (typeof (rowwise) === "undefined") {
    // check if vector of numbers
    if (type(elems) === "vector")
      return new Float64Array(elems);

    // check if 2D Array => toMatrix rowwise
    let rowwise = true;
    for (k = 0; k < elems.length; k++) {
      if (!Array.isArray(elems[k]) || elemtypes[k] === "vector") {
        rowwise = false;
        if (elemtypes[k] === "string")
          return elems; // received vector of strings => return it directly
      }
    }
  }

  if (elems.length === 0) {
    return [];
  }

  let m = 0;
  let n = 0;
  let i;
  let j;
  if (rowwise) {
    let res = [];

    for (k = 0; k < elems.length; k++) {
      switch (elemtypes[k]) {
        case "matrix":
          res.push(elems[k].val);
          m += elems[k].m;
          n = elems[k].n;
          break;

        case "vector":
          if (concatWithNumbers) {
            // return a column by concatenating vectors and numbers
            for (let l = 0; l < elems[k].length; l++) {
              res.push(elems[k][l]);
            }
            n = 1;
            m += elems[k].length;
          } else {
            // vector (auto transposed) as row in a matrix
            res.push(elems[k]);
            m += 1;
            n = elems[k].length;
          }
          break;

        case "number":
          res.push(elems[k]);
          m += 1;
          n = 1;
          break;

        case "spvector":
          return spmat(elems);

        default:
          // Array containing not only numbers...
          // probably calling mat( Array2D ) => return Array2D
          return elems;
          break;
      }
    }
    if (n === 1) {
      let M = new Float64Array(res);
      return M;
    }
    let M = new Matrix(m, n);
    let p = 0;
    for (k = 0; k < res.length; k++) {
      if (res[k].buffer) {
        M.val.set(res[k], p);
        p += res[k].length;
      } else {
        for (j = 0; j < res[k].length; j++) {
          M.val[p + j] = res[k][j];
        }
        p += res[k].length;
      }
    }
    return M;
  } else {
    // compute size
    m = size(elems[0], 1);
    for (k = 0; k < elems.length; k++) {
      if (elemtypes[k] === "matrix")
        n += elems[k].n;
      else
        n++;
      if (size(elems[k], 1) !== m)
        return "undefined";
    }

    // Build matrix
    let res = new Matrix(m, n);
    let c;
    for (i = 0; i < m; i++) {
      c = 0; // col index
      for (k = 0; k < elems.length; k++) {
        switch (elemtypes[k]) {
          case "matrix":
            for (j = 0; j < elems[k].n; j++) {
              res.val[i*n + j + c] = elems[k].val[i*elems[k].n + j];
            }
            c += elems[k].n;
            break;

          case "vector": //vector
            res.val[i*n + c] = elems[k][i];
            c++;
            break;

          case "number":
            res.val[i*n + c] = elems[k];
            c++;
            break;
          default:
            break;
        }
      }
    }

    return res;
  }
}

/// Relational Operators

export function isEqual(a, b) {
  let i;
  let j;
  let res;
  let ta = type(a);
  let tb = type(b);

  if (ta === "number" && tb !== "number")
    return isEqual(b, a);

  if (ta !== "number" && tb === "number") {
    // vector/matrix + scalar
    switch (ta) {
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (isZero(a[i] - b))
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isEqual(a.val, b), true);
        return res;
        break;
      default:
        return (a == b ? 1 : 0);
    }
  } else if (ta === tb) {

    switch (ta) {
      case "number":
        return (isZero(a - b) ? 1 : 0);
        break;
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (isZero(a[i] - b[i]))
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isEqual(a.val, b.val), true);
        return res;
        break;
      default:
        return (a == b ? 1 : 0);
    }
  } else
    return "undefined";
}

export function isNotEqual(a, b) {
  let i;
  let j;
  let res;
  let ta = type(a);
  let tb = type(b);

  if (ta === "number" && tb !== "number")
    return isNotEqual(b, a);

  if (ta !== "number" && tb === "number") {
    // vector/matrix + scalar
    switch (ta) {
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (!isZero(a[i] - b))
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isNotEqual(a.val, b), true);
        return res;
        break;
      default:
        return (a != b ? 1 : 0);
    }
  } else if (ta === tb) {

    switch (ta) {
      case "number":
        return (!isZero(a - b) ? 1 : 0);
        break;
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (!isZero(get(a, i) - get(b, i)))
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isNotEqual(a.val, b.val), true);
        return res;
        break;
      default:
        return (a != b ? 1 : 0);
    }
  } else
    return "undefined";
}

export function isGreater(a, b) {
  let i;
  let j;
  let res;
  let ta = type(a);
  let tb = type(b);

  if (ta === "number" && tb !== "number")
    return isGreater(b, a);

  if (ta !== "number" && tb === "number") {
    // vector/matrix + scalar
    switch (ta) {
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (a[i] - b > EPS)
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isGreater(a.val, b), true);
        return res;
        break;
      default:
        return (a > b ? 1 : 0);
    }
  } else if (ta === tb) {

    switch (ta) {
      case "number":
        return (a > b ? 1 : 0);
        break;
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (a[i] - b[i] > EPS)
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isGreater(a.val, b.val), true);
        return res;
        break;
      default:
        return (a > b ? 1 : 0);
    }
  } else
    return "undefined";
}

export function isGreaterOrEqual(a, b) {
  let i;
  let j;
  let res;
  let ta = type(a);
  let tb = type(b);

  if (ta === "number" && tb !== "number")
    return isGreaterOrEqual(b, a);

  if (ta !== "number" && tb === "number") {
    // vector/matrix + scalar
    switch (ta) {
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (a[i] - b > -EPS)
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isGreaterOrEqual(a.val, b), true);
        return res;
        break;
      default:
        return (a >= b ? 1 : 0);
    }
  } else if (ta === tb) {

    switch (ta) {
      case "number":
        return (a >= b);
        break;
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (a[i] - b[i] > -EPS)
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isGreaterOrEqual(a.val, b.val), true);
        return res;
        break;
      default:
        return (a >= b ? 1 : 0);
    }
  } else
    return "undefined";
}

export function isLower(a, b) {
  let i;
  let j;
  let res;
  let ta = type(a);
  let tb = type(b);

  if (ta === "number" && tb !== "number")
    return isLower(b, a);

  if (ta !== "number" && tb === "number") {
    // vector/matrix + scalar
    switch (ta) {
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (b - a[i] > EPS)
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isLower(a.val, b), true);
        return res;
        break;
      default:
        return (a < b ? 1 : 0);
    }
  } else if (ta === tb) {

    switch (ta) {
      case "number":
        return (a < b ? 1 : 0);
        break;
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (b[i] - a[i] > EPS)
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isLower(a.val, b.val), true);
        return res;
        break;
      default:
        return (a < b ? 1 : 0);
    }
  } else
    return "undefined";
}

export function isLowerOrEqual(a, b) {
  let i;
  let j;
  let res;

  let ta = type(a);
  let tb = type(b);

  if (ta === "number" && tb !== "number")
    return isLowerOrEqual(b, a);

  if (ta !== "number" && tb === "number") {
    // vector/matrix + scalar
    switch (ta) {
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (b - a[i] > -EPS)
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isLowerOrEqual(a.val, b), true);
        return res;
        break;
      default:
        return (a <= b ? 1 : 0);
    }
  } else if (ta === tb) {

    switch (ta) {
      case "number":
        return (a <= b ? 1 : 0);
        break;
      case "vector":
        res = new Float64Array(a.length);
        for (i = 0; i < a.length; i++) {
          if (b[i] - a[i] > -EPS)
            res[i] = 1;
        }
        return res;
        break;
      case "matrix":
        res = new Matrix(a.m, a.n, isLowerOrEqual(a.val, b.val), true);
        return res;
        break;
      default:
        return (a <= b ? 1 : 0);
    }
  } else
    return "undefined";
}


export function find(b) {
  // b is a boolean vector of 0 and 1.
  // return the indexes of the 1's.
  let i;
  let n = b.length;
  let res = [];
  for (i = 0; i < n; i++) {
    if (b[i] !== 0)
      res.push(i);
  }
  return res;
}

export function findmax(x) {
  // return the index of the maximum in x
  let i;

  switch (type(x)) {
    case "number":
      return 0;
    case "vector": {
      let idx = 0;
      let maxi = x[0];
      for (i = 1; i < x.length; i++) {
        if (x[i] > maxi) {
          maxi = x[i];
          idx = i;
        }
      }
      return idx;
    }
    case "spvector": {
      let maxi = x.val[0];
      let idx = x.ind[0];

      for (i = 1; i < x.val.length; i++) {
        if (x.val[i] > maxi) {
          maxi = x.val[i];
          idx = x.ind[i];
        }
      }
      if (maxi < 0 && x.val.length < x.length) {
        idx = 0;
        while (x.ind.indexOf(idx) >= 0 && idx < x.length) {
          idx++;
        }
      }
      return idx;
    }
    default:
      return "undefined";
  }

}

export const argmax = findmax;

export function findmin(x) {
  // return the index of the minimum in x
  let i;

  switch (type(x)) {
    case "number":
      return 0;
    case "vector": {
      let idx = 0;
      let mini = x[0];
      for (i = 1; i < x.length; i++) {
        if (x[i] < mini) {
          mini = x[i];
          idx = i;
        }
      }
      return idx;
    }
    case "spvector": {
      let mini = x.val[0];
      let idx = x.ind[0];

      for (i = 1; i < x.val.length; i++) {
        if (x.val[i] < mini) {
          mini = x.val[i];
          idx = x.ind[i];
        }
      }
      if (mini > 0 && x.val.length < x.length) {
        idx = 0;
        while (x.ind.indexOf(idx) >= 0 && idx < x.length) {
          idx++;
        }
      }
      return idx;
    }
    default:
      return "undefined";
  }

}

export const argmin = findmin;

/**
 * @param {Float64Array}
 * @param {boolean}
 * @param {boolean}
 * @return {Float64Array|Array}
 */
export function sort(x, decreasingOrder, returnIndexes) {
  // if returnIndexes = true : replace x with its sorted version
  // otherwise return a sorted copy without altering x

  if (typeof (decreasingOrder) === "undefined")
    decreasingOrder = false;
  if (typeof (returnIndexes) === "undefined")
    returnIndexes = false;

  let i;
  let j;
  let tmp;

  const n = x.length;
  if (returnIndexes) {
    let indexes = range(n);
    for (i = 0; i < n - 1; i++) {
      if (decreasingOrder)
        j = findmax(get(x, range(i, n))) + i;
      else
        j = findmin(get(x, range(i, n))) + i;

      if (i != j) {
        tmp = x[i];
        x[i] = x[j];
        x[j] = tmp;

        tmp = indexes[i];
        indexes[i] = indexes[j];
        indexes[j] = tmp;
      }
    }
    return indexes;
  } else {
    let xs = vectorCopy(x);
    for (i = 0; i < n - 1; i++) {
      if (decreasingOrder)
        j = findmax(get(xs, range(i, n))) + i;
      else
        j = findmin(get(xs, range(i, n))) + i;

      if (i != j) {
        tmp = xs[i];
        xs[i] = xs[j];
        xs[j] = tmp;
      }
    }
    return xs;
  }
}

/// Stats
/**
 * @param {Float64Array}
 * @return {number}
 */
export function sumVector(a) {
  let i;
  const n = a.length;
  let res = a[0];
  for (i = 1; i < n; i++) {
    res += a[i];
  }
  return res;
}

/**
 * @param {Matrix}
 * @return {number}
 */
export function sumMatrix(A) {
  return sumVector(A.val);
}

/**
 * @param {Matrix}
 * @return {Matrix}
 */
export function sumMatrixRows(A) {
  let i;
  let j;
  const m = A.m;
  const n = A.n;
  let res = new Float64Array(n);
  let r = 0;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      res[j] += A.val[r + j];
    }
    r += n;
  }
  return new Matrix(1, n, res, true); // return row vector
}

/**
 * @param {Matrix}
 * @return {Float64Array}
 */
export function sumMatrixCols(A) {
  const m = A.m;
  let res = new Float64Array(m);
  let r = 0;
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < A.n; j++) {
      res[i] += A.val[r + j];
    }
    r += A.n;
  }
  return res;
}

export function sum(A, sumalongdimension) {

  switch (type(A)) {
    case "vector":
      if (arguments.length === 1 || sumalongdimension === 1) {
        return sumVector(A);
      } else {
        return vectorCopy(A);
      }
      break;
    case "spvector":
      if (arguments.length === 1 || sumalongdimension === 1)
        return sumVector(A.val);
      else
        return A.copy();
      break;

    case "matrix":
      if (arguments.length === 1) {
        return sumMatrix(A);
      } else if (sumalongdimension === 1) {
        return sumMatrixRows(A);
      } else if (sumalongdimension === 2) {
        return sumMatrixCols(A);
      } else
        return undefined;
      break;
    case "spmatrix":
      if (arguments.length === 1) {
        return sumVector(A.val);
      } else if (sumalongdimension === 1) {
        return sumspMatrixRows(A);
      } else if (sumalongdimension === 2) {
        return sumspMatrixCols(A);
      } else
        return undefined;
      break;
    default:
      return A;
      break;
  }
}

/**
 * @param {Float64Array}
 * @return {number}
 */
export function prodVector(a) {
  let i;
  const n = a.length;
  let res = a[0];
  for (i = 1; i < n; i++) {
    res *= a[i];
  }
  return res;
}

/**
 * @param {Matrix}
 * @return {number}
 */
export function prodMatrix(A) {
  return prodVector(A.val);
}

/**
 * @param {Matrix}
 * @return {Matrix}
 */
export function prodMatrixRows(A) {
  let i;
  let j;
  const m = A.m;
  const n = A.n;
  let res = new Float64Array(A.row(0));
  let r = n;
  for (i = 1; i < m; i++) {
    for (j = 0; j < n; j++) {
      res[j] *= A.val[r + j];
    }
    r += A.n;
  }
  return new Matrix(1, n, res, true); // return row vector
}

/**
 * @param {Matrix}
 * @return {Float64Array}
 */
export function prodMatrixCols(A) {
  const m = A.m;
  let res = new Float64Array(m);
  let r = 0;
  for (let i = 0; i < m; i++) {
    res[i] = A.val[r];
    for (let j = 1; j < A.n; j++) {
      res[i] *= A.val[r + j];
    }
    r += A.n;
  }
  return res;
}

export function prod(A, prodalongdimension) {

  switch (type(A)) {
    case "vector":
      if (arguments.length === 1 || prodalongdimension === 1)
        return prodVector(A);
      else
        return vectorCopy(A);
    case "spvector":
      if (arguments.length === 1 || prodalongdimension === 1) {
        if (A.val.length < A.length)
          return 0;
        else
          return prodVector(A.val);
      } else
        return A.copy();
    case "matrix":
      if (arguments.length === 1) {
        return prodMatrix(A);
      } else if (prodalongdimension === 1) {
        return prodMatrixRows(A);
      } else if (prodalongdimension === 2) {
        return prodMatrixCols(A);
      } else
        return undefined;
    case "spmatrix":
      if (arguments.length === 1) {
        if (A.val.length < A.m*A.n)
          return 0;
        else
          return prodVector(A.val);
      } else if (prodalongdimension === 1) {
        return prodspMatrixRows(A);
      } else if (prodalongdimension === 2) {
        return prodspMatrixCols(A);
      } else
        return undefined;
    default:
      return A;
  }
}

export function mean(A, sumalongdimension) {

  switch (type(A)) {
    case "vector":
      if (arguments.length === 1 || sumalongdimension === 1) {
        return sumVector(A)/A.length;
      } else {
        return vectorCopy(A);
      }
    case "spvector":
      if (arguments.length === 1 || sumalongdimension === 1)
        return sumVector(A.val)/A.length;
      else
        return A.copy();

    case "matrix":
      if (arguments.length === 1) {
        return sumMatrix(A)/(A.m*A.n);
      } else if (sumalongdimension === 1) {
        return mulScalarMatrix(1/A.m, sumMatrixRows(A));
      } else if (sumalongdimension === 2) {
        return mulScalarVector(1/A.n, sumMatrixCols(A));
      } else
        return undefined;
    case "spmatrix":
      if (arguments.length === 1) {
        return sumVector(A.val)/(A.m*A.n);
      } else if (sumalongdimension === 1) {
        return mulScalarMatrix(1/A.m, sumspMatrixRows(A));
      } else if (sumalongdimension === 2) {
        return mulScalarVector(1/A.n, sumspMatrixCols(A));
      } else
        return undefined;
    default:
      return A;
  }
}

export function variance(A, alongdimension) {
  let meanA;

  // variance = sum(A^2)/n - mean(A)^2
  if (arguments.length > 1)
    meanA = mean(A, alongdimension);
  else
    meanA = mean(A);

  switch (type(A)) {
    case "number":
      return 0;
    case "vector":
      if (arguments.length === 1 || alongdimension === 1) {
        let res = (dot(A, A)/A.length) - meanA*meanA;
        return res;
      } else {
        return zeros(A.length);
      }
    case "spvector":
      if (arguments.length === 1 || alongdimension === 1) {
        let res = (dot(A.val, A.val)/A.length) - meanA*meanA;
        return res;
      } else
        return zeros(A.length);

    case "matrix":
    case "spmatrix":
      if (typeof (alongdimension) === "undefined") {
        let res = (sum(entrywisemul(A, A))/(A.m*A.n)) - meanA*meanA;
        return res;
      } else if (alongdimension === 1) {
        // let of columns
        let res = sub(entrywisediv(sum(entrywisemul(A, A), 1), A.length), entrywisemul(meanA, meanA));
        return res;
      } else if (alongdimension === 2) {
        // sum all columns, result is column vector
        let res = sub(entrywisediv(sum(entrywisemul(A, A), 2), A.n), entrywisemul(meanA, meanA));
        return res;
      } else
        return undefined;
    default:
      return undefined;
  }
}

export function std(A, alongdimension) {
  if (arguments.length > 1)
    return sqrt(variance(A, alongdimension));
  else
    return sqrt(variance(A));
}

/**
 * Covariance matrix C = X'*X ./ X.m
 * @param {Matrix|Float64Array|spVector}
 * @return {Matrix|number}
 */
export function cov(X) {
  switch (type(X)) {
    case "number":
      return 0;
    case "vector": {
      let mu = mean(X);
      return (dot(X, X)/X.length - mu*mu);
    }
    case "spvector": {
      let mu = mean(X);
      return (dot(X.val, X.val)/X.length - mu*mu);
    }
    case "matrix": {
      let mu = mean(X, 1).row(0);
      return divMatrixScalar(xtx(subMatrices(X, outerprod(ones(X.m), mu))), X.m);
    }
    case "spmatrix": {
      let mu = mean(X, 1).row(0);
      return divMatrixScalar(xtx(subspMatrixMatrix(X, outerprod(ones(X.m), mu))), X.m);
    }
    default:
      return undefined;
  }
}

/**
 * Compute X'*X
 * @param {Matrix}
 * @return {Matrix}
 */
export function xtx(X) {
  const N = X.m;
  const d = X.n;

  let C = new Matrix(d, d);
  for (let i = 0; i < N; i++) {
    let xi = X.row(i);
    for (let k = 0; k < d; k++) {
      let xik = xi[k];
      for (let j = k; j < d; j++) {
        C.val[k*d + j] += xik*xi[j];
      }
    }
  }
  // Symmetric lower triangular part:
  for (let k = 0; k < d; k++) {
    let kd = k*d;
    for (let j = k; j < d; j++) {
      C.val[j*d + k] = C.val[kd + j];
    }
  }
  return C;
}

export function norm(A, sumalongdimension) {
  // l2-norm (Euclidean norm) of vectors or Frobenius norm of matrix
  let i;
  let j;
  switch (type(A)) {
    case "number":
      return Math.abs(A);
      break;
    case "vector":
      if (arguments.length === 1 || sumalongdimension === 1) {
        return Math.sqrt(dot(A, A));
      } else
        return abs(A);
      break;
    case "spvector":
      if (arguments.length === 1 || sumalongdimension === 1) {
        return Math.sqrt(dot(A.val, A.val));
      } else
        return abs(A);
      break;
    case "matrix":
      if (arguments.length === 1) {
        return Math.sqrt(dot(A.val, A.val));
      } else if (sumalongdimension === 1) {
        // norm of columns, result is row vector
        const n = A.n;
        let res = zeros(1, n);
        let r = 0;
        for (i = 0; i < A.m; i++) {
          for (j = 0; j < n; j++) {
            res.val[j] += A.val[r + j]*A.val[r + j];
          }
          r += n;
        }
        for (j = 0; j < n; j++) {
          res.val[j] = Math.sqrt(res.val[j]);
        }
        return res;
      } else if (sumalongdimension === 2) {
        // norm of rows, result is column vector
        let res = zeros(A.m);
        let r = 0;
        for (i = 0; i < A.m; i++) {
          for (j = 0; j < A.n; j++) {
            res[i] += A.val[r + j]*A.val[r + j];
          }
          r += A.n;
          res[i] = Math.sqrt(res[i]);
        }

        return res;
      } else
        return "undefined";
      break;
    case "spmatrix":
      if (arguments.length === 1) {
        return Math.sqrt(dot(A.val, A.val));
      } else if (sumalongdimension === 1 && !A.rowmajor) {
        // norm of columns, result is row vector
        const nn = A.n;
        let res = zeros(1, nn);
        for (j = 0; j < nn; j++) {
          let s = A.cols[j];
          let e = A.cols[j + 1];
          for (let k = s; k < e; k++) {
            res.val[j] += A.val[k]*A.val[k];
          }
          res.val[j] = Math.sqrt(res.val[j]);
        }
        return res;
      } else if (sumalongdimension === 2 && A.rowmajor) {
        // norm of rows, result is column vector
        let res = zeros(A.m);
        for (i = 0; i < A.m; i++) {
          let s = A.rows[i];
          let e = A.rows[i + 1];
          for (let k = s; k < e; k++) {
            res[i] += A.val[k]*A.val[k];
          }
          res[i] = Math.sqrt(res[i]);
        }

        return res;
      } else
        return "undefined";
      break;
    default:
      return "undefined";
  }
}

export function norm1(A, sumalongdimension) {
  // l1-norm of vectors and matrices
  if (arguments.length === 1)
    return sum(abs(A));
  else
    return sum(abs(A), sumalongdimension);
}

export function norminf(A, sumalongdimension) {
  // linf-norm of vectors and max-norm of matrices
  if (arguments.length === 1)
    return max(abs(A));
  else
    return max(abs(A), sumalongdimension);
}

export function normp(A, p, sumalongdimension) {
  // lp-norm of vectors and matrices
  if (arguments.length === 2)
    return Math.pow(sum(pow(abs(A), p)), 1/p);
  else
    return pow(sum(pow(abs(A), p), sumalongdimension), 1/p);
}

export function normnuc(A) {
  // nuclear norm
  switch (type(A)) {
    case "matrix":
      return sumVector(svd(A));
      break;
    case "spmatrix":
      return sumVector(svd(fullMatrix(A)));
      break;
    case "number":
      return A;
      break;
    case "vector":
    case "spvector":
      return 1;
      break;
    default:
      return undefined;
      break;
  }
}

export function norm0(A, sumalongdimension, epsilonarg) {
  // l0-pseudo-norm of vectors and matrices
  // if epsilon > 0, consider values < epsilon as 0

  let epsilon = EPS;
  if (arguments.length === 3)
    epsilon = epsilonarg;

  let i;
  let j;
  switch (type(A)) {
    case "number":
      return (Math.abs(A) > epsilon);
      break;
    case "vector":
      if (arguments.length === 1 || sumalongdimension === 1) {
        return norm0Vector(A, epsilon);
      } else
        return isGreater(abs(a), epsilon);
      break;
    case "spvector":
      if (arguments.length === 1 || sumalongdimension === 1) {
        return norm0Vector(A.val, epsilon);
      } else
        return isGreater(abs(a), epsilon);
      break;
    case "matrix":
      if (arguments.length === 1) {
        return norm0Vector(A.val, epsilon);
      } else if (sumalongdimension === 1) {
        // norm of columns, result is row vector
        let res = zeros(1, A.n);
        for (i = 0; i < A.m; i++) {
          for (j = 0; j < A.n; j++) {
            if (Math.abs(A[i*A.n + j]) > epsilon)
              res.val[j]++;
          }
        }
        return res;
      } else if (sumalongdimension === 2) {
        // norm of rows, result is column vector
        let res = zeros(A.m);
        for (i = 0; i < A.m; i++) {
          for (j = 0; j < A.n; j++) {
            if (Math.abs(A[i*A.n + j]) > epsilon)
              res[i]++;
          }
        }
        return res;
      } else
        return undefined;
      break;
    case "spmatrix":
      if (arguments.length === 1) {
        return norm0Vector(A.val, epsilon);
      } else if (sumalongdimension === 1) {
        // norm of columns, result is row vector
        let res = zeros(1, A.n);
        if (A.rowmajor) {
          for (let k = 0; k < A.val.length; k++) {
            if (Math.abs(A.val[k]) > epsilon)
              res.val[A.cols[k]]++;
          }
        } else {
          for (let i = 0; i < A.n; i++) {
            res.val[i] = norm0Vector(A.col(i).val, epsilon);
          }
        }
        return res;
      } else if (sumalongdimension === 2) {
        // norm of rows, result is column vector
        let res = zeros(A.m);
        if (A.rowmajor) {
          for (let i = 0; i < A.m; i++) {
            res[i] = norm0Vector(A.row(i).val, epsilon);
          }
        } else {
          for (let k = 0; k < A.val.length; k++) {
            if (Math.abs(A.val[k]) > epsilon)
              res[A.rows[k]]++;
          }
        }
        return res;
      } else
        return undefined;
      break;
    default:
      return undefined;
  }
}

/**
 * @param {Float64Array}
 * @param {number}
 * @return {number}
 */
export function norm0Vector(x, epsilon) {
  const n = x.length;
  let res = 0;
  for (let i = 0; i < n; i++) {
    if (Math.abs(x[i]) > epsilon)
      res++;
  }
  return res;
}

///////////////////////////////////////////:
// Linear systems of equations
///////////////////////////////////////

export function solve(A, b) {
  /* Solve the linear system Ax = b	*/

  let tA = type(A);

  if (tA === "vector" || tA === "spvector" || (tA === "matrix" && A.m === 1)) {
    // One-dimensional least squares problem:
    let AtA = mul(transpose(A), A);
    let Atb = mul(transpose(A), b);
    return Atb/AtA;
  }

  if (tA === "spmatrix") {
    /*if ( A.m === A.n )
			return spsolvecg(A, b); // assume A is positive definite
		else*/
    return spcgnr(A, b);
  }

  if (type(b) === "vector") {
    if (A.m === A.n)
      return solveGaussianElimination(A, b);
    else
      return solveWithQRcolumnpivoting(A, b);
  } else
    return solveWithQRcolumnpivotingMultipleRHS(A, b); // b is a matrix
}

/**
 * Solve the linear system Ax = b given the Cholesky factor L of A
 * @param {Matrix}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function cholsolve(L, b) {
  let z = forwardsubstitution(L, b);
  let x = backsubstitution(transposeMatrix(L), z);
  return x;
}

/**
 * @param {Matrix}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function solveWithQRfactorization(A, b) {
  const m = A.length;
  const n = A.n;
  let QRfact = qr(A);
  let R = QRfact.R;
  let beta = QRfact.beta;

  let btmp = vectorCopy(b);
  let j;
  let i;
  let k;
  let v;

  let smallb;

  for (j = 0; j < n - 1; j++) {
    v = get(R, range(j, m), j); // get Householder vectors
    v[0] = 1;
    // b(j:m) = (I - beta v v^T ) * b(j:m)
    smallb = get(btmp, range(j, m));
    set(btmp, range(j, m), sub(smallb, mul(beta[j]*mul(v, smallb), v)));
  }
  // last iteration only if m>n
  if (m > n) {
    j = n - 1;

    v = get(R, range(j, m), j); // get Householder vectors
    v[0] = 1;
    // b(j:m) = (I - beta v v^T ) * b(j:m)
    smallb = get(btmp, range(j, m));
    set(btmp, range(j, m), sub(smallb, mul(beta[j]*mul(v, smallb), v)));

  }

  // Solve R x = b with backsubstitution (R is upper triangular, well it is not really here because we use the lower part to store the vectors v):
  return backsubstitution(R, get(btmp, range(n)));


//	return backsubstitution ( get ( R, range(n), range(n) ) , rows ( btmp, range(1,n)) );
//	we can spare the get and copy of R : backsubstitution will only use this part anyway
}

/**
 * @param {Matrix}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function backsubstitution(U, b) {
  // backsubstitution to solve a linear system U x = b with upper triangular U

  const n = b.length;
  let j = n - 1;
  let x = zeros(n);

  if (!isZero(U.val[j*n + j]))
    x[j] = b[j]/U.val[j*n + j];

  j = n - 2;
  if (!isZero(U.val[j*n + j]))
    x[j] = (b[j] - U.val[j*n + n - 1]*x[n - 1])/U.val[j*n + j];

  for (j = n - 3; j >= 0; j--) {
    if (!isZero(U.val[j*n + j]))
      x[j] = (b[j] - dot(U.row(j).subarray(j + 1, n), x.subarray(j + 1, n)))/U.val[j*n + j];
  }

  // solution
  return x;
}

/**
 * @param {Matrix}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function forwardsubstitution(L, b) {
  // forward substitution to solve a linear system L x = b with lower triangular L

  const n = b.length;
  let j;
  let x = zeros(n);

  if (!isZero(L.val[0]))
    x[0] = b[0]/L.val[0];

  if (!isZero(L.val[n + 1]))
    x[1] = (b[1] - L.val[n]*x[0])/L.val[n + 1];

  for (j = 2; j < n; j++) {
    if (!isZero(L.val[j*n + j]))
      x[j] = (b[j] - dot(L.row(j).subarray(0, j), x.subarray(0, j)))/L.val[j*n + j];
  }

  // solution
  return x;
}

/**
 * @param {Matrix}
 * @param {Float64Array}
 * @return {Float64Array}
 */
export function solveWithQRcolumnpivoting(A, b) {

  let m;
  let n;
  let R;
  let V;
  let beta;
  let r;
  let piv;
  if (type(A) === "matrix") {
    // Compute the QR factorization
    m = A.m;
    n = A.n;
    let QRfact = qr(A);
    R = QRfact.R;
    V = QRfact.V;
    beta = QRfact.beta;
    r = QRfact.rank;
    piv = QRfact.piv;
  } else {
    // we get the QR factorization in A
    R = A.R;
    r = A.rank;
    V = A.V;
    beta = A.beta;
    piv = A.piv;
    m = R.m;
    n = R.n;
  }

  let btmp = vectorCopy(b);
  let j;
  let i;
  let k;

  let smallb;
  // b = Q' * b
  for (j = 0; j < r; j++) {

    // b(j:m) = (I - beta v v^T ) * b(j:m)
    smallb = get(btmp, range(j, m));

    set(btmp, range(j, m), sub(smallb, mul(beta[j]*mul(V[j], smallb), V[j])));
  }
  // Solve R x = b with backsubstitution
  let x = zeros(n);

  if (r > 1) {
    set(x, range(0, r), backsubstitution(R, get(btmp, range(r))));
    // note: if m < n, backsubstitution only uses n columns of R.
  } else {
    x[0] = btmp[0]/R.val[0];
  }

  // and apply permutations
  for (j = r - 1; j >= 0; j--) {
    if (piv[j] !== j) {
      let tmp = x[j];
      x[j] = x[piv[j]];
      x[piv[j]] = tmp;
    }
  }
  return x;

}

/**
 * @param {Matrix}
 * @param {Matrix}
 * @return {Matrix}
 */
export function solveWithQRcolumnpivotingMultipleRHS(A, B) {

  let m;
  let n;
  let R;
  let V;
  let beta;
  let r;
  let piv;
  if (type(A) === "matrix") {
    // Compute the QR factorization
    m = A.m;
    n = A.n;
    let QRfact = qr(A);
    R = QRfact.R;
    V = QRfact.V;
    beta = QRfact.beta;
    r = QRfact.rank;
    piv = QRfact.piv;
  } else {
    // we get the QR factorization in A
    R = A.R;
    r = A.rank;
    V = A.V;
    beta = A.beta;
    piv = A.piv;
    m = R.m;
    n = R.n;
  }

  let btmp = matrixCopy(B);
  let j;
  let i;
  let k;

  let smallb;
  // B = Q' * B
  for (j = 0; j < r; j++) {

    // b(j:m) = (I - beta v v^T ) * b(j:m)
    smallb = get(btmp, range(j, m), []);

    set(btmp, range(j, m), [], sub(smallb, mul(mul(beta[j], V[j]), mul(transpose(V[j]), smallb))));
  }
  // Solve R X = B with backsubstitution
  let X = zeros(n, m);

  if (r > 1) {
    for (j = 0; j < m; j++) {
      set(X, range(0, r), j, backsubstitution(R, get(btmp, range(r), j)));
    }
    // note: if m < n, backsubstitution only uses n columns of R.
  } else {
    set(X, 0, [], entrywisediv(get(btmp, 0, []), R.val[0]));
  }

  // and apply permutations
  for (j = r - 1; j >= 0; j--) {
    if (piv[j] !== j) {
      swaprows(X, j, piv[j]);
    }
  }
  return X;

}

export function solveGaussianElimination(Aorig, borig) {

  // Solve square linear system Ax = b with Gaussian elimination

  let i;
  let j;
  let k;

  let A = matrixCopy(Aorig).toArrayOfFloat64Array(); // useful to quickly switch rows
  let b = vectorCopy(borig);

  const m = Aorig.m;
  const n = Aorig.n;
  if (m !== n)
    return undefined;

  // Set to zero small values... ??

  for (k = 0; k < m; k++) {

    // Find imax = argmax_i=k...m |A_i,k|
    let imax = k;
    let Aimaxk = Math.abs(A[imax][k]);
    for (i = k + 1; i < m; i++) {
      let Aik = Math.abs(A[i][k]);
      if (Aik > Aimaxk) {
        imax = i;
        Aimaxk = Aik;
      }
    }
    if (isZero(Aimaxk)) {
      console.log("** Warning in solve(A,b), A is square but singular, switching from Gaussian elimination to QR method.");
      return solveWithQRcolumnpivoting(Aorig, borig);
    }

    if (imax !== k) {
      // Permute the rows
      let a = A[k];
      A[k] = A[imax];
      A[imax] = a;
      let tmpb = b[k];
      b[k] = b[imax];
      b[imax] = tmpb;
    }
    let Ak = A[k];

    // Normalize row k
    let Akk = Ak[k];
    b[k] /= Akk;

    //Ak[k] = 1; // not used afterwards
    for (j = k + 1; j < n; j++) {
      Ak[j] /= Akk;
    }

    if (Math.abs(Akk) < 1e-8) {
      console.log("** Warning in solveGaussianElimination: " + Akk + " " + k + ":" + m);
    }

    // Substract the kth row from others to get 0s in kth column
    let Aik;
    let bk = b[k];
    for (i = 0; i < m; i++) {
      if (i !== k) {
        let Ai = A[i];
        Aik = Ai[k];
        for (j = k + 1; j < n; j++) { // Aij = 0  with j < k and Aik = 0 after this operation but is never used
          Ai[j] -= Aik*Ak[j];
        }
        b[i] -= Aik*bk;
      }
    }
  }

  // Solution:
  return b;
}

export function inv(M) {
  if (typeof (M) === "number")
    return 1/M;

  // inverse matrix with Gaussian elimination

  let i;
  let j;
  let k;
  const m = M.length;
  const n = M.n;
  if (m !== n)
    return "undefined";

  // Make extended linear system:
  let A = matrixCopy(M);
  let B = eye(n);

  for (k = 0; k < m; k++) {
    let kn = k*n;

    // Find imax = argmax_i=k...m |A_i,k|
    let imax = k;
    let Aimaxk = Math.abs(A.val[imax*n + k]);
    for (i = k + 1; i < m; i++) {
      if (Math.abs(A.val[i*n + k]) > Aimaxk) {
        imax = i;
        Aimaxk = Math.abs(A.val[i*n + k]);
      }
    }
    if (Math.abs(Aimaxk) < 1e-12) {
      return "singular";
    }

    if (imax !== k) {
      // Permute the rows
      swaprows(A, k, imax);
      swaprows(B, k, imax);
    }

    // Normalize row k
    let Akk = A.val[kn + k];
    for (j = 0; j < n; j++) {
      A.val[kn + j] /= Akk;
      B.val[kn + j] /= Akk;
    }

    if (Math.abs(Akk) < 1e-8)
      console.log("!! Warning in inv(): " + Akk + " " + k + ":" + m);

    // Substract the kth row from others to get 0s in kth column
    let Aik;
    for (i = 0; i < m; i++) {
      if (i !== k) {
        let ri = i*n;
        Aik = A.val[ri + k];
        if (!isZero(Aik)) {
          for (j = 0; j < n; j++) {
            A.val[ri + j] -= Aik*A.val[kn + j];
            B.val[ri + j] -= Aik*B.val[kn + j];
          }
        }
      }
    }
  }

  // Solution:
  return B;
}

export function chol(A) {
  // Compute the Cholesky factorization A = L L^T with L lower triangular
  // for a positive definite and symmetric A
  // returns L or undefined if A is not positive definite
  const n = A.m;
  if (A.n !== n) {
    error("Cannot compute the cholesky factorization: the matrix is not square.");
    return undefined;
  }
  const n2 = n*n;
  const Aval = A.val;
  let L = new Float64Array(n2);

  let i, j;
  // first column = A(:,0) / sqrt(L(0,0)
  let sqrtLjj = Math.sqrt(Aval[0]);
  for (i = 0; i < n2; i += n) { 	// i = i*n = ptr to row i
    L[i] = Aval[i]/sqrtLjj;
  }
  // other colums
  j = 1;
  let jn = n;
  while (j < n && !isNaN(sqrtLjj)) {
    for (i = jn; i < n2; i += n) {	// i = i*n
      let Lij = Aval[i + j];
      for (let k = 0; k < j; k++) {
        Lij -= L[jn + k]*L[i + k];
      }
      if (i === jn)
        sqrtLjj = Math.sqrt(Lij);

      L[i + j] = Lij/sqrtLjj;
    }
    j++;
    jn += n;
  }
  if (isNaN(sqrtLjj))
    return undefined; // not positive definite
  else
    return new Matrix(n, n, L, true);
}

export function ldlsymmetricpivoting(Aorig) {
  // LDL factorization for symmetric matrices
  let A = matrixCopy(Aorig);
  let n = A.length;
  if (A.m !== n) {
    error("Error in ldl(): the matrix is not square.");
    return undefined;
  }
  let k;
  let piv = zeros(n);
  let alpha;
  let v;

  for (k = 0; k < n - 1; k++) {

    piv[k] = findmax(get(diag(A), range(k, n)));
    swaprows(A, k, piv[k]);
    swapcols(A, k, piv[k]);
    alpha = A.val[k*n + k];
    v = getCols(A, [k]).subarray(k + 1, n);

    for (let i = k + 1; i < n; i++) {
      A.val[i*n + k] /= alpha;
    }

    set(A, range(k + 1, n), range(k + 1, n), sub(get(A, range(k + 1, n), range(k + 1, n)), outerprod(v, v, 1/alpha)));

  }

  // Make it lower triangular
  for (let j = 0; j < n - 1; j++) {
    for (let k = j + 1; k < n; k++) {
      A.val[j*n + k] = 0;
    }
  }
  return {L: A, piv: piv};
}

/**
 * @param {Float64Array}
 * @return {{v: Float64Array, beta: number}}
 */
export function house(x) {
  // Compute Houselholder vector v such that
  // P = (I - beta v v') is orthogonal and Px = ||x|| e_1

  const n = x.length;
  let i;
  let mu;
  let beta;
  let v = zeros(n);
  let v0;
  let sigma;

  let x0 = x[0];
  let xx = dot(x, x);

  // sigma = x(2:n)^T x(2:n)
  sigma = xx - x0*x0;

  if (isZero(sigma)) {
    // x(2:n) is zero =>  v=[1,0...0], beta = 0
    beta = 0;
    v[0] = 1;
  } else {
    mu = Math.sqrt(xx); // norm(x) ; //Math.sqrt( x0*x0 + sigma );
    if (x0 < EPS) {
      v0 = x0 - mu;
    } else {
      v0 = -sigma/(x0 + mu);
    }

    beta = 2*v0*v0/(sigma + v0*v0);

    // v = [v0,x(2:n)] / v0
    v[0] = 1;
    for (i = 1; i < n; i++) {
      v[i] = x[i]/v0;
    }
  }

  return {"v": v, "beta": beta};
}

/**
 * @param {Matrix}
 * @return {{Q: (Matrix|undefined), R: Matrix, beta: Float64Array}
 */
export function qroriginal(A, compute_Q) {
  // QR factorization based on Householder reflections WITHOUT column pivoting
  // A with m rows and n cols; m >= n

  // test with A = [[12,-51,4],[6,167,-68],[-4,24,-41]]
  // then R = [ [14 -21 -14 ], [ -3, 175, -70], [2, -0.75, 35]]

  let m = A.length;
  let n = A.n;
  if (n > m)
    return "QR factorization unavailable for n > m.";

  let i;
  let j;
  let k;
  let householder;
  let R = matrixCopy(A);
  let beta = zeros(n);
  let outer;
  let smallR;
  let Q;
  let V = []; // store householder vectors


  for (j = 0; j < n - 1; j++) {
    householder = house(get(R, range(j, m), j));
    // R(j:m,j:n) = ( I - beta v v' ) * R(j:m,j:n) = R - (beta v) (v'R)
    smallR = get(R, range(j, m), range(j, n));
    set(R, range(j, m), range(j, n), subMatrices(smallR, outerprodVectors(householder.v, mulMatrixVector(transposeMatrix(smallR), householder.v), householder.beta)));

    V[j] = householder.v;
    beta[j] = householder.beta;

  }
  // Last iteration only if m > n: if m=n, (I - beta v v' ) = 1 => R(n,n) is unchanged
  if (m > n) {
    j = n - 1;
    smallR = get(R, range(j, m), j)
    householder = house(smallR);
    // R(j:m,n) = ( I - beta v v' ) * R(j:m, n) = R(j:m,n) - (beta v) (v'R(j:m,n) ) = Rn - ( beta *(v' * Rn) )* v
    set(R, range(j, m), n - 1, subVectors(smallR, mulScalarVector(dot(householder.v, smallR)*householder.beta, householder.v)));

    V[j] = vectorCopy(householder.v);
    beta[j] = householder.beta;

  }

  if (compute_Q) {
    let r;
    if (typeof (compute_Q) === "number") {
      // compute only first r columns of Q
      r = compute_Q;
      Q = eye(m, r);
    } else {
      Q = eye(m);
      r = m;
    }
    let smallQ;
    let nmax = n - 1;
    if (m <= n)
      nmax = n - 2;
    if (nmax >= r)
      nmax = r - 1;

    for (j = nmax; j >= 0; j--) {
      smallQ = get(Q, range(j, m), range(j, r));

      if (r > 1) {
        if (j === r - 1)
          set(Q, range(j, m), [j], subVectors(smallQ, mulScalarVector(dot(smallQ, V[j])*beta[j], V[j])));
        else
          set(Q, range(j, m), range(j, r), sub(smallQ, outerprod(V[j], mul(transpose(smallQ), V[j]), beta[j])));
      } else
        Q = subVectors(smallQ, mulScalarVector(dot(smallQ, V[j])*beta[j], V[j]));
    }
  }

  return {"Q": Q, "R": R, "beta": beta};
}

/**
 * @param {Matrix}
 * @return {{Q: (Matrix|undefined), R: Matrix, V: Array, beta: Float64Array, piv: Float64Array, rank: number}
 */
export function qr(A, compute_Q) {
  // QR factorization with column pivoting AP = QR based on Householder reflections
  // A with m rows and n cols; m >= n (well, it also works with m < n)
  // piv = vector of permutations : P = P_rank with P_j = identity with swaprows ( j, piv(j) )

  // Implemented with R transposed for faster computations on rows instead of columns

  /* TEST
	A  = [[12,-51,4],[6,167,-68],[-4,24,-41]]
	QR = qr(A)
	QR.R


	*/
  const m = A.m;
  const n = A.n;

  /*
	if ( n > m)
		return "QR factorization unavailable for n > m.";
	*/

  let i;
  let j;

  let householder;
  let R = transpose(A);// transposed for faster implementation
  let Q;

  let V = []; // store householder vectors in this list (not a matrix)
  let beta = zeros(n);
  let piv = zeros(n);

  let smallR;

  let r = -1; // rank estimate -1

  let normA = norm(A);
  let normR22 = normA;
  let Rij;

  const TOL = 1e-5;
  let TOLnormR22square = TOL*normA;
  TOLnormR22square *= TOLnormR22square;

  let tau = 0;
  let k = 0;
  let c = zeros(n);
  for (j = 0; j < n; j++) {
    let Rj = R.val.subarray(j*R.n, j*R.n + R.n);
    c[j] = dot(Rj, Rj);
    if (c[j] > tau) {
      tau = c[j];
      k = j;
    }
  }

  let updateR = function (r, v, beta) {
    // set ( R, range(r,n), range(r,m) , subMatrices (  smallR , outerprodVectors( mulMatrixVector( smallR, householder.v), householder.v,  householder.beta ) ) ) ;
    // most of the time is spent here...
    let i, j, l;
    let m_r = m - r;
    for (i = r; i < n; i++) {
      let smallRiv = 0;
      let Ri = i*m + r; // =  i * R.n + r
      let Rval = R.val.subarray(Ri, Ri + m_r);
      for (l = 0; l < m_r; l++) {
        smallRiv += Rval[l]*v[l];
      }	//smallRiv += R.val[Ri + l] * v[l];
      smallRiv *= beta;
      for (j = 0; j < m_r; j++) {
        Rval[j] -= smallRiv*v[j]; // R.val[Ri + j] -= smallRiv * v[j];
      }
    }
  };

  // Update c
  let updateC = function (r) {
    let j;
    for (j = r + 1; j < n; j++) {
      let Rjr = R.val[j*m + r];
      c[j] -= Rjr*Rjr;
    }

    // tau, k = max ( c[r+1 : n] )
    k = r + 1;
    tau = c[r + 1];
    for (j = r + 2; j < n; j++) {
      if (c[j] > tau) {
        tau = c[j];
        k = j;
      }
    }
  };

  // Compute norm of residuals
  let computeNormR22 = function (r) {
    //normR22 = norm(get ( R, range(r+1,n), range(r+1,m), ) );
    let normR22 = 0;
    let i = r + 1;
    let ri = i*m;
    let j;
    while (i < n && normR22 <= TOLnormR22square) {
      for (j = r + 1; j < m; j++) {
        let Rij = R.val[ri + j];
        normR22 += Rij*Rij;
      }
      i++;
      ri += m;
    }
    return normR22;
  }


  while (tau > EPS && r < n - 1 && normR22 > TOLnormR22square) {

    r++;

    piv[r] = k;
    swaprows(R, r, k);
    c[k] = c[r];
    c[r] = tau;

    if (r < m - 1) {
      householder = house(R.val.subarray(r*R.n + r, r*R.n + m)); // house only reads vec so subarray is ok
    } else {
      householder.v = [1];
      householder.beta = 0;
      //smallR = R[m-1][m-1];
    }

    if (r < n - 1) {
      // smallR is a matrix
      updateR(r, householder.v, householder.beta);
    } else {
      // smallR is a row vector (or a number if m=n):
      if (r < m - 1) {
        updateR(r, householder.v, householder.beta);
        /*
				let r_to_m = range(r,m);
				smallR = get(R, r, r_to_m);
				set ( R, r , r_to_m, sub (  smallR , transpose(mul( householder.beta * mul( smallR, householder.v) ,householder.v  ) )) ) ;*/
      } else {
        //let smallRnumber = R.val[(m-1)*R.n + m-1]; // beta is zero, so no update
        //set ( R, r , r, sub (  smallRnumber , transpose(mul( householder.beta * mul( smallRnumber, householder.v) ,householder.v  ) )) ) ;
      }
    }

    // Store householder vectors and beta
    V[r] = vectorCopy(householder.v);
    beta[r] = householder.beta;

    if (r < n - 1) {
      // Update c
      updateC(r);

      // stopping criterion for rank estimation
      if (r < m - 1)
        normR22 = computeNormR22(r);
      else
        normR22 = 0;
    }
  }

  if (compute_Q) {
    Q = eye(m);
    let smallQ;
    let nmax = r;
    if (m > r + 1)
      nmax = r - 1;
    for (j = nmax; j >= 0; j--) {
      if (j === m - 1) {
        Q.val[j*m + j] -= beta[j]*V[j][0]*V[j][0]*Q.val[j*m + j];
      } else {
        let j_to_m = range(j, m);
        smallQ = get(Q, j_to_m, j_to_m);// matrix
        set(Q, j_to_m, j_to_m, subMatrices(smallQ, outerprodVectors(V[j], mulMatrixVector(transposeMatrix(smallQ), V[j]), beta[j])));
      }
    }
  }

  return {"Q": Q, "R": transpose(R), "V": V, "beta": beta, "piv": piv, "rank": r + 1};
}

export function qrRnotTransposed(A, compute_Q) {
  // QR factorization with column pivoting AP = QR based on Householder reflections
  // A with m rows and n cols; m >= n (well, it also works with m < n)
  // piv = vector of permutations : P = P_rank with P_j = identity with swaprows ( j, piv(j) )

  // original implementation working on columns

  /* TEST
	A  = [[12,-51,4],[6,167,-68],[-4,24,-41]]
	QR = qr(A)
	QR.R


	*/
  let m = A.m;
  let n = A.n;

  /*
	if ( n > m)
		return "QR factorization unavailable for n > m.";
	*/

  let i;
  let j;

  let householder;
  let R = matrixCopy(A);
  let Q;

  let V = []; // store householder vectors in this list (not a matrix)
  let beta = zeros(n);
  let piv = zeros(n);

  let smallR;

  let r = -1; // rank estimate -1

  let normA = norm(A);
  let normR22 = normA;

  let TOL = 1e-6;

  let tau = 0;
  let k = 0;
  let c = zeros(n);
  for (j = 0; j < n; j++) {
    let Aj = getCols(A, [j]);
    c[j] = dot(Aj, Aj);
    if (c[j] > tau) {
      tau = c[j];
      k = j;
    }
  }

  while (tau > EPS && r < n - 1 && normR22 > TOL*normA) {

    r++;

    piv[r] = k;
    swapcols(R, r, k);
    c[k] = c[r];
    c[r] = tau;

    if (r < m - 1) {
      householder = house(get(R, range(r, m), r));
      smallR = get(R, range(r, m), range(r, n));
    } else {
      householder.v = [1];
      householder.beta = 0;
      smallR = R[m - 1][m - 1];
    }

    if (r < n - 1) {
      // smallR is a matrix
      set(R, range(r, m), range(r, n), subMatrices(smallR, outerprodVectors(householder.v, mulMatrixVector(transposeMatrix(smallR), householder.v), householder.beta)));
    } else {
      // smallR is a vector (or a number if m=n):
      set(R, range(r, m), r, sub(smallR, mul(householder.beta*mul(smallR, householder.v), householder.v)));
    }

    // Store householder vectors and beta
    if (m > r + 1)
      V[r] = vectorCopy(householder.v);
    beta[r] = householder.beta;

    if (r < n - 1) {
      // Update c
      for (j = r + 1; j < n; j++) {
        c[j] -= R[r][j]*R[r][j];
      }

      // tau, k = max ( c[r+1 : n] )
      k = r + 1;
      tau = c[r + 1];
      for (j = r + 2; j < n; j++) {
        if (c[j] > tau) {
          tau = c[j];
          k = j;
        }
      }

      // stopping criterion for rank estimation
      if (r < m - 1) {
        //normR22 = norm(get ( R, range(r+1,m),range(r+1,n) ) );
        normR22 = 0;
        for (i = r + 1; i < m; i++) {
          for (j = r + 1; j < n; j++) {
            let Rij = R[i][j];
            normR22 += Rij*Rij;
          }
        }
        normR22 = Math.sqrt(normR22);
      } else
        normR22 = 0;
    }
  }

  if (compute_Q) {
    Q = eye(m);
    let smallQ;
    let nmax = r;
    if (m > r + 1)
      nmax = r - 1;
    for (j = nmax; j >= 0; j--) {
      if (j === m - 1) {
        Q.val[j*m + j] -= beta[j]*V[j][0]*V[j][0]*Q.val[j*m + j];
      } else {
        smallQ = get(Q, range(j, m), range(j, m));
        set(Q, range(j, m), range(j, m), subMatrices(smallQ, outerprodVectors(V[j], mulMatrixVector(transposeMatrix(smallQ), V[j]), beta[j])));
      }
    }

  }

  return {"Q": Q, "R": R, "V": V, "beta": beta, "piv": piv, "rank": r + 1};
}

/** Conjugate gradient method for solving the symmetyric positive definite system Ax = b
 * @param{{Matrix|spMatrix}}
 * @param{Float64Array}
 * @return{Float64Array}
 */
export function solvecg(A, b) {
  if (A.type === "spmatrix")
    return spsolvecg(A, b);
  else
    return solvecgdense(A, b);
}

/** Conjugate gradient method for solving the symmetyric positive definite system Ax = b
 * @param{Matrix}
 * @param{Float64Array}
 * @return{Float64Array}
 */
export function solvecgdense(A, b) {
  /*
TEST
A = randn(2000,1000)
x = randn(1000)
b = A*x + 0.01*randn(2000)
tic()
xx = solve(A,b)
t1 = toc()
ee = norm(A*xx - b)
tic()
xh=solvecg(A'*A, A'*b)
t2 = toc()
e = norm(A*xh - b)
*/

  const n = A.n;
  const m = A.m;

  let x = randn(n); //vectorCopy(x0);
  let r = subVectors(b, mulMatrixVector(A, x));
  let rhoc = dot(r, r);
  const TOL = 1e-8;
  let delta2 = TOL*norm(b);
  delta2 *= delta2;

  // first iteration:
  let p = vectorCopy(r);
  let w = mulMatrixVector(A, p);
  let mu = rhoc/dot(p, w);
  saxpy(mu, p, x);
  saxpy(-mu, w, r);
  let rho_ = rhoc;
  rhoc = dot(r, r);

  let k = 1;

  let updateP = function (tau, r) {
    for (let i = 0; i < m; i++) {
      p[i] = r[i] + tau*p[i];
    }
  }

  while (rhoc > delta2 && k < n) {
    updateP(rhoc/rho_, r);
    w = mulMatrixVector(A, p);
    mu = rhoc/dot(p, w);
    saxpy(mu, p, x);
    saxpy(-mu, w, r);
    rho_ = rhoc;
    rhoc = dot(r, r);
    k++;
  }
  return x;
}

/** Conjugate gradient normal equation residual method for solving the rectangular system Ax = b
 * @param{{Matrix|spMatrix}}
 * @param{Float64Array}
 * @return{Float64Array}
 */
export function cgnr(A, b) {
  if (A.type === "spmatrix")
    return spcgnr(A, b);
  else
    return cgnrdense(A, b);
}

/** Conjugate gradient normal equation residual method for solving the rectangular system Ax = b
 * @param{Matrix}
 * @param{Float64Array}
 * @return{Float64Array}
 */
export function cgnrdense(A, b) {
  /*
TEST
A = randn(2000,1000)
x = randn(1000)
b = A*x + 0.01*randn(2000)
tic()
xx = solve(A,b)
t1 = toc()
ee = norm(A*xx - b)
tic()
xh=cgnr(A, b)
t2 = toc()
e = norm(A*xh - b)
*/

  const n = A.n;
  const m = A.m;

  let x = randn(n); // vectorCopy(x0);
  let At = transposeMatrix(A);
  let r = subVectors(b, mulMatrixVector(A, x));
  const TOL = 1e-8;
  let delta2 = TOL*norm(b);
  delta2 *= delta2;

  // first iteration:
  let z = mulMatrixVector(At, r);
  let rhoc = dot(z, z);
  let p = vectorCopy(z);
  let w = mulMatrixVector(A, p);
  let mu = rhoc/dot(w, w);
  saxpy(mu, p, x);
  saxpy(-mu, w, r);
  z = mulMatrixVector(At, r);
  let rho_ = rhoc;
  rhoc = dot(z, z);

  let k = 1;

  let updateP = function (tau, z) {
    for (let i = 0; i < m; i++) {
      p[i] = z[i] + tau*p[i];
    }
  }

  while (rhoc > delta2 && k < n) {
    updateP(rhoc/rho_, z);
    w = mulMatrixVector(A, p);
    mu = rhoc/dot(w, w);
    saxpy(mu, p, x);
    saxpy(-mu, w, r);
    z = mulMatrixVector(At, r);
    rho_ = rhoc;
    rhoc = dot(z, z);
    k++;
  }
  return x;
}

/** Lanczos algorithm
 * @param{Matrix}
 */
export function lanczos(A, q1) {

  const maxIters = 300;
  const TOL = EPS*norm(A);
  const n = A.n;
  let i;
  let k = 0;
  let w = vectorCopy(q1);
  let v = mulMatrixVector(A, w);
  let alpha = dot(w, v);

  saxpy(-alpha, w, v);

  let beta = norm(b);

  while (beta > TOL && k < maxIters) {

    for (i = 0; i < n; i++) {
      let t = w[i];
      w[i] = v[i]/beta;
      v[i] = -beta/t;
    }

    let Aw = mulMatrixVector(A, w);

    for (i = 0; i < n; i++) {
      v[i] += Aw[i];
    }

    alpha = dot(w, v);
    saxpy(-alpha, w, v);
    beta = norm(v);
    k++;
  }
}

/**
 * @param{Matrix}
 * @param{boolean}
 * @return{Matrix}
 */
export function tridiagonalize(A, returnQ) {
  // A : a square and symmetric  matrix
  // T = Q A Q' , where T is tridiagonal and Q = (H1 ... Hn-2)' is the product of Householder transformations.
  // if returnQ, then T overwrites A
  let k;
  const n = A.length;
  let T;
  let Q;
  let Pk;
  if (returnQ) {
    T = A;
    Q = eye(n);
    let beta = [];
    let V = [];
  } else
    T = matrixCopy(A);
  let p;
  let w;
  let vwT;
  let normTkp1k;
  let householder;

  for (k = 0; k < n - 2; k++) {
    let Tkp1k = get(T, range(k + 1, n), k);
    let Tkp1kp1 = get(T, range(k + 1, n), range(k + 1, n));

    householder = house(Tkp1k);
    p = mulScalarVector(householder.beta, mulMatrixVector(Tkp1kp1, householder.v));
    w = subVectors(p, mulScalarVector(0.5*householder.beta*dot(p, householder.v), householder.v));

    /*
		T[k+1][k] = norm ( Tkp1k );
		T[k][k+1] = T[k+1][k];
		*/
    // make T really tridiagonal: the above does not modify the other entries to set them to 0
    normTkp1k = zeros(n - k - 1);
    normTkp1k[0] = norm(Tkp1k);
    set(T, k, range(k + 1, n), normTkp1k);
    set(T, range(k + 1, n), k, normTkp1k);

    vwT = outerprodVectors(householder.v, w);
    set(T, range(k + 1, n), range(k + 1, n), subMatrices(subMatrices(Tkp1kp1, vwT), transpose(vwT)));

    if (returnQ) {
      V[k] = householder.v;
      beta[k] = householder.beta;
    }
  }
  if (returnQ) {
    let updateQ = function (j, v, b) {
      // Q = Q - b* v (Q'v)'
      //smallQ =  get(Q, range(j,n), range(j,n) );// matrix
      //set ( Q, range(j,n), range(j,n) , subMatrices (  smallQ , outerprodVectors(  V[k], mulMatrixVector( transposeMatrix(smallQ), V[k]), beta[k] ) ) );
      let i, k;
      let Qtv = zeros(n - j);
      let n_j = n - j;
      for (i = 0; i < n_j; i++) {
        let Qi = (i + j)*n + j;
        for (k = 0; k < n_j; k++) {
          Qtv[k] += v[i]*Q.val[Qi + k];
        }
      }
      for (i = 0; i < n_j; i++) {
        let Qi = (i + j)*n + j;
        let betavk = b*v[i];
        for (k = 0; k < n_j; k++) {
          Q.val[Qi + k] -= betavk*Qtv[k];
        }
      }
    };

    // Backaccumulation of Q
    for (k = n - 3; k >= 0; k--) {
      updateQ(k + 1, V[k], beta[k]);
    }
    return Q;
  } else
    return T;
}

export function givens(a, b, Gi, Gk, n) {
  // compute a Givens rotation:
  let c;
  let s;
  let tau;
  let G;

  // Compute c and s
  if (b === 0) {
    c = 1;
    s = 0;
  } else {
    if (Math.abs(b) > Math.abs(a)) {
      tau = -a/b;
      s = 1/Math.sqrt(1 + tau*tau);
      c = s*tau;
    } else {
      tau = -b/a;
      c = 1/Math.sqrt(1 + tau*tau);
      s = c*tau;
    }
  }

  if (arguments.length === 5) {
    // Build Givens matrix G from c and s:
    G = eye(n);
    G.val[Gi*n + Gi] = c;
    G.val[Gi*n + Gk] = s;
    G.val[Gk*n + Gi] = -s;
    G.val[Gk*n + Gk] = c;
    return G;
  } else {
    return [c, s];
  }

}

/**
 * @param {number}
 * @param {number}
 * @param {number}
 * @param {number}
 * @param {Matrix}
 */
export function premulGivens(c, s, i, k, A) {
  // apply a Givens rotation to A : A([i,k],:) = G' * A([i,k],:)
  //  with G = givens (a,b,i,k) and [c,s]=givens(a,b)
  // NOTE: this modifies A

  const n = A.n;
  let j;
  const ri = i*n;
  const rk = k*n;
  let t1;
  let t2;
  for (j = 0; j < n; j++) {
    t1 = A.val[ri + j];
    t2 = A.val[rk + j];
    A.val[ri + j] = c*t1 - s*t2;
    A.val[rk + j] = s*t1 + c*t2;
  }
}

/**
 * @param {number}
 * @param {number}
 * @param {number}
 * @param {number}
 * @param {Matrix}
 */
export function postmulGivens(c, s, i, k, A) {
  // apply a Givens rotation to A : A(:, [i,k]) =  A(:, [i,k]) * G
  //  with G = givens (a,b,i,k) and [c,s]=givens(a,b)
  // NOTE: this modifies A

  const m = A.length;
  let j;
  let t1;
  let t2;
  let rj = 0;
  for (j = 0; j < m; j++) {
    t1 = A.val[rj + i];
    t2 = A.val[rj + k];
    A.val[rj + i] = c*t1 - s*t2;
    A.val[rj + k] = s*t1 + c*t2;
    rj += A.n;
  }
}

export function implicitSymQRWilkinsonShift(T, computeZ) {
  // compute T = Z' T Z
  // if computeZ:  return {T,cs} such that T = Z' T Z  with Z = G1.G2...
  // and givens matrices Gk of parameters cs[k]

  const n = T.length;
  const rn2 = n*(n - 2);
  const rn1 = n*(n - 1);

  const d = (T.val[rn2 + n - 2] - T.val[rn1 + n - 1])/2;
  const t2 = T.val[rn1 + n - 2]*T.val[rn1 + n - 2];
  const mu = T.val[rn1 + n - 1] - t2/(d + Math.sign(d)*Math.sqrt(d*d + t2));
  let x = T.val[0] - mu; // T[0][0]
  let z = T.val[n];		// T[1][0]
  let cs, csArray;

  if (computeZ) {
    csArray = new Array(n - 1);
  }
  //let Z = eye(n);

  let k;
  for (k = 0; k < n - 1; k++) {
    /*
		G = givens(x,z, k, k+1, n);
		T = mul(transpose(G), mul(T, G) ); // can do this much faster
		if ( computeZ ) {
			Z = mul(Z, G );
		}
		*/
    cs = givens(x, z);
    postmulGivens(cs[0], cs[1], k, k + 1, T);
    premulGivens(cs[0], cs[1], k, k + 1, T);
    if (computeZ)
      csArray[k] = [cs[0], cs[1]];
    //postmulGivens(cs[0], cs[1], k, k+1, Z);

    if (k < n - 2) {
      let r = n*(k + 1) + k;
      x = T.val[r];
      z = T.val[r + n]; // [k+2][k];
    }
  }
  if (computeZ) {
    return {"T": T, "cs": csArray};
//		return {"T": T, "Z": Z} ;
  } else
    return T;
}

export function eig(A, computeEigenvectors) {
  // Eigendecomposition of a symmetric matrix A (QR algorithm)

  let Q;
  let D;
  if (computeEigenvectors) {
    D = matrixCopy(A);
    Q = tridiagonalize(D, true);
  } else {
    D = tridiagonalize(A);
  }

  let q;
  let p;
  const n = A.length;
  let i;

  const TOL = 1e-12; //10 * EPS;

  do {
    for (i = 0; i < n - 1; i++) {
      if (Math.abs(D.val[i*n + i + 1]) < TOL*(Math.abs(D.val[i*n + i]) + Math.abs(D.val[(i + 1)*n + i + 1]))) {
        D.val[i*n + i + 1] = 0;
        D.val[(i + 1)*n + i] = 0;
      }
    }

    // find largest q such that D[n-p-q:n][n-p-q:n] is diagonal:
    if (!isZero(D.val[(n - 1)*n + n - 2]) || !isZero(D.val[(n - 2)*n + n - 1]))
      q = 0;
    else {
      q = 1;
      while (q < n - 1 && isZero(D.val[(n - q - 1)*n + n - q - 2]) && isZero(D.val[(n - q - 2)*n + n - q - 1])) {
        q++;
      }
      if (q >= n - 1)
        q = n;
    }

    // find smallest p such that D[p:q][p:q] is unreduced ( without zeros on subdiagonal?)
    p = -1;
    let zerosOnSubdiagonal;
    do {
      p++;
      zerosOnSubdiagonal = false;
      let k = p;

      while (k < n - q - 1 && zerosOnSubdiagonal === false) {
        if (isZero(D.val[(k + 1)*n + k]))
          zerosOnSubdiagonal = true;
        k++;
      }
    } while (zerosOnSubdiagonal && p + q < n);

    // Apply implicit QR iteration
    if (q < n) {

      if (computeEigenvectors) {
        let res = implicitSymQRWilkinsonShift(get(D, range(p, n - q), range(p, n - q)), true);
        set(D, range(p, n - q), range(p, n - q), res.T);
        for (let kk = 0; kk < n - q - p - 1; kk++) {
          postmulGivens(res.cs[kk][0], res.cs[kk][1], p + kk, p + kk + 1, Q);
        }
        //Z = eye(n);
        //set(Z, range(p,n-q), range(p,n-q), DZ22.Z );
        // Q = mulMatrixMatrix ( Q, Z );

      } else {
        set(D, range(p, n - q), range(p, n - q), implicitSymQRWilkinsonShift(get(D, range(p, n - q), range(p, n - q)), false));
      }
    }

  } while (q < n) ;

  if (computeEigenvectors) {
    return {"V": diag(D), "U": Q};
  } else
    return diag(D);
}

export function eigs(A, r, smallest) {
  // Compute r largest or smallest eigenvalues and eigenvectors
  if (typeof (r) === "undefined")
    r = 1;
  if (typeof (smallest) === "undefined" || smallest === false || smallest != "smallest") {
    if (r === 1)
      return eig_powerIteration(A);
    else
      return eig_orthogonalIteration(A, r);
  } else {
    // look for smallest eigenvalues
    if (r === 1)
      return eig_inverseIteration(A, 0);
    else
      return eig_bisect(A, r);
    //return eig_inverseOrthogonalIteration ( A , r) ;
  }
}

export function eig_powerIteration(A, u0) {
// Compute the largest eigenvalue and eigenvector with the power method
  const maxIters = 1000;
  let k;
  const n = A.length;

  // init with a random u or an initial guess u0
  let u;
  if (typeof (u0) === "undefined")
    u = randn(n);
  else
    u = u0;
  u = mulScalarVector(1/norm(u), u);
  let lambda = 1;
  for (k = 0; k < maxIters; k++) {
    // Apply the iteration : u = Au / norm(Au)
    u = mulMatrixVector(A, u);
    lambda = norm(u);
    u = mulScalarVector(1/lambda, u);
  }
  return {"v": lambda, "u": u};
}

export function eig_orthogonalIteration(A, r) {

  if (r === 1)
    return eig_powerIteration(A);

// Compute the r largest eigenvalue and eigenvector with the power method (orthogonal iteration)
  const maxIters = 1000;
  let k;
  const n = A.length;

  // init with a random Q
  let Q = randn(n, r);
  let normQ = norm(Q, 1);
  Q = entrywisediv(Q, mul(ones(n), normQ));
  let QR;
  let Z;

  const TOL = 1e-11;
  let V;

  for (k = 0; k < maxIters; k++) {

    // Z = AQ
    Z = mulMatrixMatrix(A, Q);
    if (Math.floor(k/50) === k/50) {
      // convergence test
      V = mulMatrixMatrix(transpose(Q), Z);

      if (norm(subMatrices(Z, mulMatrixMatrix(Q, diag(diag(V))))) < TOL)
        break;
    }

    // QR = Z	// XXX maybe not do this at every iteration...
    Q = qroriginal(Z, r).Q;

  }

  V = mulMatrixMatrix(transpose(Q), mulMatrixMatrix(A, Q));

  return {"V": diag(V), "U": Q};
}

export function eig_inverseIteration(A, lambda) {
  // Compute an eigenvalue-eigenvector pair from an approximate eigenvalue with the inverse iteration
  let perturbation = 0.0001*lambda;

  let maxIters;

  if (typeof (maxIters) === "undefined")
    maxIters = 100;

  let k;
  const n = A.length;

  // apply power iteration with (A - lambda I)^-1 instead of A
  let A_lambdaI = sub(A, mul(lambda + perturbation, eye(n)));
  let QR = qr(A_lambdaI); // and precompute QR factorization

  while (QR.rank < n) { // check if not singular
    perturbation *= 10;
    A_lambdaI = sub(A, mul(lambda + perturbation, eye(n)));
    QR = qr(A_lambdaI); // and precompute QR factorization
    //console.log(perturbation);
  }

  // init
  let u = sub(mul(2, rand(n)), 1); //ones(n); //
  u = mulScalarVector(1/norm(u), u);
  let v;
  let r;
  let norminfA = norminf(A);
  k = 0;
  do {
    // u =  solve(A_lambdaI , u) ;

    u = solveWithQRcolumnpivoting(QR, u); // QR factorization precomputed

    v = norm(u);
    u = entrywisediv(u, v);

    r = mulMatrixVector(A_lambdaI, u);

    k++;
  } while (k < maxIters && maxVector(absVector(r)) < 1e-10*norminfA); // && Math.abs(v * perturbation - 1 ) < EPS );
  return u;

}

export function eigenvector(A, lambda) {
  return eig_inverseIteration(A, lambda, 2);
}

export function eig_inverseOrthogonalIteration(A, r) {

  if (r === 1)
    return eig_inverseIteration(A);

// Compute the r smallest eigenvalue and eigenvectors with the inverse power method
// (orthogonal iteration)
  const maxIters = 1000;
  let k;
  const n = A.length;
  let QR = qr(A); // precompute QR factorization

  // init with a random Q
  let Q = randn(n, r);
  let normQ = norm(Q, 1);
  Q = entrywisediv(Q, mul(ones(n), normQ));
  let Z;

  const TOL = 1e-11;
  let V;

  for (k = 0; k < maxIters; k++) {

    // Z = A^-1 Q
    Z = solveWithQRcolumnpivotingMultipleRHS(QR, Q);

    if (Math.floor(k/50) === k/50) {
      // convergence test
      V = mulMatrixMatrix(transpose(Q), Z);

      if (norm(subMatrices(Z, mulMatrixMatrix(Q, V))) < TOL)
        break;
    }

    // QR = Z	// XXX maybe not do this at every iteration...
    Q = qroriginal(Z, r).Q;

  }

  V = mulMatrixMatrix(transpose(Q), mulMatrixMatrix(A, Q));

  return {"V": diag(V), "U": Q, "iters": k};
}


export function eig_bisect(A, K) {
// find K smallest eigenvalues

  /*
TEST
//Symmetric eigenvalue decomposition
X = rand(5,5)
A = X*X'
v = eig(A)
eig_bisect(A,3)
*/

  let x, y, z;

  // Tridiagonalize A
  let T = tridiagonalize(A);
  const n = T.n;
  let a = diag(T);
  let b = zeros(n);
  let i;
  for (i = 0; i < n - 1; i++) {
    b[i] = T.val[i*n + i + 1];
  }

  // Initialize [y,z] with Gershgorin disk theorem
  let y0 = a[0] - b[0];
  let z0 = a[0] + b[0];
  for (let i = 1; i < n; i++) {
    let yi = a[i] - b[i] - b[i - 1];
    let zi = a[i] + b[i] + b[i - 1];
    if (yi < y0)
      y0 = yi;
    if (zi > z0)
      z0 = zi;
  }

  /*
	// polynomial evaluation and counting sign changes (original method)
	let polya = function (x,a,b,n) {
		let pr_2 = 1;
		let pr_1 = a[0] - x;
		let pr;
		let signchanges = 0;
		if (  pr_1 < EPS )
			signchanges = 1;

		let r;
		for ( r = 1; r < n ; r++) {
			pr = (a[r] - x) * pr_1 - b[r-1] * b[r-1] * pr_2;

			if ( Math.abs(pr) < EPS || (pr > 0 &&  pr_1 < 0 ) || (pr < 0) && (pr_1 > 0) )
				signchanges ++;

			pr_2 = pr_1;
			pr_1 = pr;
		}
		return signchanges;
	};
	*/

  // ratio of polynomials evaluation and counting sign changes
  // (modification discussed in Barth et al., 1967 for better stability due to pr ~ 0 in the above)
  let polyq = function (x, a, b, n) {
    let qi_1 = a[0] - x;
    let qi;
    let signchanges = 0;
    if (qi_1 < EPS)
      signchanges = 1;

    let i;
    for (i = 1; i < n; i++) {
      qi = (a[i] - x) - b[i - 1]*b[i - 1]/qi_1;

      if (qi < EPS)
        signchanges++;

      if (Math.abs(qi) < EPS)
        qi_1 = EPS;
      else
        qi_1 = qi;
    }
    return signchanges;
  };


  // Start bisection
  const TOL = 1e-10;
  let lambda = zeros(K);
  let xu = entrywisemul(z0, ones(K)); // upper bounds on lambdas
  y = y0;
  let n_lowerthan_x;// nb of eigenvalues lower than x
  for (let k = 1; k <= K; k++) {
    // k is the number of desired eigenvalues in this sweep

    z = xu[k - 1];
    //y=y; from previous sweep

    // find the (n-k+1)th eigenvalue
    while (Math.abs(z - y) > TOL*(Math.abs(y) + Math.abs(z))) {
      x = (y + z)/2;
      n_lowerthan_x = polyq(x, a, b, n);

      if (n_lowerthan_x >= k)
        z = x; // enough eigenvalues below x, decrease upper bound to x
      else
        y = x; // not enough ev below x, increase lower bound to x

      // update boudns on other lambdas
      for (let j = k + 1; j <= K; j++) {
        if (n_lowerthan_x >= j)
          xu[j - 1] = x;
      }

    }
    lambda[k - 1] = (y + z)/2;
  }
  //return lambda;

  // Compute eigenvectors: XXX can be faster by using inverse iteration on the tridiagonal matrix
  //						 with faster system solving

  let u = eigenvector(A, lambda[0]);
  let U = mat([u], false);

  for (let k = 1; k < K; k++) {
    // deal with too close eigenvalues
    let perturbtol = 10*Math.max(EPS, Math.abs(EPS*lambda[k - 1]));
    if (lambda[k] < lambda[k - 1] + perturbtol)
      lambda[k] = lambda[k - 1] + perturbtol;

    u = eigenvector(A, lambda[k]);
    U = mat([U, u], false);
    U = qroriginal(U, U.n).Q; // orthogonalize
  }


  return {U: U, V: lambda};
}


export function bidiagonalize(A, computeU, thinU, computeV) {
  // B = U' A V , where B is upper bidiagonal

  let j;
  const m = A.length;
  const n = A.n;
  let B;
  B = matrixCopy(A);

  let householder;

  if (computeU) {
    if (thinU) {
      let U = eye(m, n);
      let nU = n;
    } else {
      let U = eye(m);
      let nU = m;
    }
  }
  if (computeV) {
    let V = eye(n);
  }


  let updateB1 = function (j, v, beta) {
    // B = B - (beta v) ( v'* B) = B-outer(beta v, B'*v)
    //Bjmjn = get ( B, range(j,m), range(j, n));
    //set ( B, range(j,m), range(j,n), sub ( Bjmjn , outerprod ( householder.v, mul(transpose(Bjmjn), householder.v), householder.beta) ) );

    let i, k;
    let Btv = zeros(n - j);
    let n_j = n - j;
    let m_j = m - j;
    for (i = 0; i < m_j; i++) {
      let Bi = (i + j)*n + j;
      for (k = 0; k < n_j; k++) {
        Btv[k] += v[i]*B.val[Bi + k];
      }
    }
    for (i = 0; i < m_j; i++) {
      let betavk = beta*v[i];
      let Bi = (i + j)*n + j;
      for (k = 0; k < n_j; k++) {
        B.val[Bi + k] -= betavk*Btv[k];
      }
    }
  };
  let updateB2 = function (j, v, beta) {
    // B = B - beta (Bv) v' (with B = B_j:m, j+1:n)

    //Bjmjn = get ( B, range(j,m), range(j+1, n));
    //set ( B, range(j,m), range(j+1,n) , sub( Bjmjn, outerprod( mul(Bjmjn, householder.v), householder.v, householder.beta) ) );
    let i, k;
    let n_j_1 = n - j - 1;
    for (i = j; i < m; i++) {
      let Bi = i*n + j + 1;
      let Bv = 0;
      for (k = 0; k < n_j_1; k++) {
        Bv += B.val[Bi + k]*v[k];
      }
      let betaBvk = beta*Bv;
      for (k = 0; k < n_j_1; k++) {
        B.val[Bi + k] -= betaBvk*v[k];
      }
    }
  };

  if (computeV) {
    let updateV = function (j, v, beta) {
      //smallV = get ( V, range(0,n), range(j+1, n));
      //set ( V, range(0,n), range(j+1,n) , sub( smallV, outerprod( mul(smallV, householder.v), householder.v, householder.beta) ) );
      let i, k;
      let n_j_1 = n - j - 1;
      for (i = 0; i < n; i++) {
        let Vi = i*n + j + 1;
        let Vv = 0;
        for (k = 0; k < n_j_1; k++) {
          Vv += V.val[Vi + k]*v[k];
        }
        let betaVvk = beta*Vv;
        for (k = 0; k < n_j_1; k++) {
          V.val[Vi + k] -= betaVvk*v[k];
        }
      }
    };
  }
  if (computeU) {
    let hv = new Array(n);// Householder vectors and betas
    let hb = new Array(n);
  }

  for (j = 0; j < n; j++) {

    if (j < m - 1) {
      householder = house(get(B, range(j, m), j));

      updateB1(j, householder.v, householder.beta);

      if (computeU) {
        hv[j] = vectorCopy(householder.v);
        hb[j] = householder.beta;
        //	updateU(j, householder.v, householder.beta);
      }
    }

    if (j < n - 2) {
      householder = house(B.row(j).subarray(j + 1, n));

      updateB2(j, householder.v, householder.beta);

      if (computeV) {
        updateV(j, householder.v, householder.beta);

      }
    }
  }
  if (computeU) {
    // Back accumulation of U (works with less than m columns)
    // Un_1 = (I-beta v v')Un = Un - beta v (v' Un)

    /*for (j=n-1;j>=0; j--) {
			if (j<m-1){
				smallU = get(U,range(j,m),[]);
				set(U,range(j,m),[], sub(smallU, mul(bv[j],mul(hv[j], mul(transpose(hv[j]) , smallU)))));
			}
		}*/
    let updateU = function (j, v, beta) {
      let i, k;
      let vtU = zeros(nU);
      for (i = j; i < m; i++) {
        let Ui = i*nU;
        let i_j = i - j;
        for (k = 0; k < nU; k++) {
          vtU[k] += v[i_j]*U.val[Ui + k];
        }
      }
      for (i = j; i < m; i++) {
        let betavk = beta*v[i - j];
        let Ui = i*nU;
        for (k = 0; k < nU; k++) {
          U.val[Ui + k] -= betavk*vtU[k];
        }
      }
    };
    let nj = Math.min(n - 1, m - 2);
    for (j = nj; j >= 0; j--) {
      updateU(j, hv[j], hb[j]);
    }
  }

  if (computeU && computeV) {
    return {"U": U, "V": V, "B": B};
  } else if (computeV)
    return {"V": V, "B": B};
  else if (computeU)
    return {"U": U, "B": B};
  else
    return B;
}


export function GolubKahanSVDstep(B, i, j, m, n, computeUV) {
  // Apply GolubKahanSVDstep to B(i:i+m, j:j+n)
  // Note: working on Utrans
  if (type(B) !== "matrix")
    return B;

  if (n < 2)
    return B;

  const rn2 = (i + n - 2)*B.n + j;
  const dm = B.val[rn2 + n - 2];
  const fm = B.val[rn2 + n - 1];
  let fm_1;
  if (n > 2)
    fm_1 = B.val[(i + n - 3)*B.n + j + n - 2];
  else
    fm_1 = 0;

  const dn = B.val[(i + n - 1)*B.n + j + n - 1];

  const d = (dm*dm + fm_1*fm_1 - dn*dn - fm*fm)/2;
  const t2 = dm*fm*dm*fm;
  const mu = dn*dn + fm*fm - t2/(d + Math.sign(d)*Math.sqrt(d*d + t2));

  let k;

  //let B0 = getCols ( B, [0]);
  //let B1 = getCols ( B, [1]) ;
  //let y = mul( B0, B0 ) - mu;
  //let z =  mul( B0, B1 );
  let y = -mu;
  let z = 0.0;
  let r0 = i*B.n + j;
  for (k = 0; k < n; k++) {
    y += B.val[r0]*B.val[r0];
    z += B.val[r0]*B.val[r0 + 1];
    r0 += B.n;
  }


  let G;
  let cs;

  let postmulgivens = function (c, s, k1, k2) {
    // apply a Givens rotation to a subset of rows of B : B(i:i+m, [k1,k2]) =  B(i:i+m, [k1,k2]) * G
    let jj;
    let t1;
    let t2;
    let rj = i*B.n + j;
    for (jj = 0; jj < m; jj++) {
      t1 = B.val[rj + k1];
      t2 = B.val[rj + k2];
      B.val[rj + k1] = c*t1 - s*t2;
      B.val[rj + k2] = s*t1 + c*t2;
      rj += B.n;
    }
  }
  let premulgivens = function (c, s, k1, k2) {
    // apply a Givens rotation to a subset of cols of B : B([k1,k2],j:j+n) = G' * B([k1,k2],j:j+n)
    let jj;
    const ri = (i + k1)*B.n + j;
    const rk = (i + k2)*B.n + j;
    let t1;
    let t2;
    for (jj = 0; jj < n; jj++) {
      t1 = B.val[ri + jj];
      t2 = B.val[rk + jj];
      B.val[ri + jj] = c*t1 - s*t2;
      B.val[rk + jj] = s*t1 + c*t2;
    }
  }

  if (computeUV) {
    //let U = eye(m);
    //let V = eye(n);
    let csU = new Array(n - 1);
    let csV = new Array(n - 1);
  }

  for (k = 0; k < n - 1; k++) {
    cs = givens(y, z);
    postmulgivens(cs[0], cs[1], k, k + 1);

    if (computeUV) {
      csV[k] = [cs[0], cs[1]];
      //	postmulGivens(cs[0],cs[1], k, k+1, V);
    }


    y = B.val[(i + k)*B.n + j + k];
    z = B.val[(i + k + 1)*B.n + j + k];

    cs = givens(y, z);
    premulgivens(cs[0], cs[1], k, k + 1);

    if (computeUV) {
      csU[k] = [cs[0], cs[1]];
      //premulGivens(cs[0],cs[1], k, k+1, U);
    }

    if (k < n - 2) {
      y = B.val[(i + k)*B.n + j + k + 1];
      z = B.val[(i + k)*B.n + j + k + 2];
    }

  }

  if (computeUV)
    return {csU: csU, csV: csV};
}

export function svd(A, computeUV) {
  /* TEST:
A=[ [-149,-50,-154],[537,180,546],[-27,-9,-25]]
s=svd(A)
should return [ 817.7597, 2.4750, 0.0030]
*/

  if (type(A) === "vector" || (type(A) === "matrix" && A.n === 1)) {
    return {"U": matrixCopy(A), "S": ones(1, 1), "V": ones(1, 1), "s": [1]};
  }
  if (A.m === 1) {
    return {"U": ones(1, 1), "S": ones(1, 1), "V": transpose(A), "s": [1]};
  }


  let i;
  let m = A.length;
  let n = A.n;


  let Atransposed = false;
  if (n > m) {
    Atransposed = true;
    let At = transposeMatrix(A);
    n = m;
    m = At.length;
  }

  let U, V, Vt, B;

  let computeU = false;
  let computeV = false;
  let thinU = false;
  if (typeof (computeUV) !== "undefined" && computeUV !== false) {

    if (computeUV === "full") {
      computeU = true;
      computeV = true;
      thinU = false;
    } else if (computeUV === true || computeUV === "thin") {
      computeU = true;
      computeV = true;
      thinU = true;
    } else if (typeof (computeUV) === "string") {
      if (computeUV.indexOf("U") >= 0)
        computeU = true;
      if (computeUV.indexOf("V") >= 0)
        computeV = true;
      if (computeUV.indexOf("thin") >= 0)
        thinU = true;
    }
    let UBV;
    if (Atransposed) {
      let tmp = computeU;
      computeU = computeV;
      computeV = tmp;
      UBV = bidiagonalize(At, computeU, thinU, computeV);
    } else
      UBV = bidiagonalize(A, computeU, thinU, computeV);

    if (computeU) {
      U = transpose(UBV.U);//Utrans
    } else
      U = undefined;

    if (computeV) {
      V = UBV.V;
      Vt = transposeMatrix(V);
    } else
      V = undefined;

    B = UBV.B;
  } else {
    if (Atransposed)
      B = bidiagonalize(At, false, false, false);
    else
      B = bidiagonalize(matrixCopy(A), false, false, false);
  }

  let B22;
  let U22;
  let V22;
  let cs;

  let q;
  let p;
  let k;

  const TOL = 1e-11;
  let iter = 0;
  do {

    for (i = 0; i < n - 1; i++) {
      if (Math.abs(B.val[i*B.n + i + 1]) < TOL*(Math.abs(B.val[i*B.n + i]) + Math.abs(B.val[(i + 1)*B.n + i + 1]))) {
        B.val[i*B.n + i + 1] = 0;
      }
    }

    // find largest q such that B[n-q+1:n][n-q+1:n] is diagonal (in matlab notation):
    q = 0;
    while (q < n && Math.abs(B.val[(n - q - 1)*B.n + n - q - 2]) < TOL && Math.abs(B.val[(n - q - 2)*B.n + n - q - 1]) < TOL) {
      q++;
    }
    if (q === n - 1)
      q = n;

    // find smallest p such that B[p+1:n-q][p+1:n-q] has no zeros on superdiag (in matlab notation):
    p = 0;	// size of B11 = first index of B22 in our notation
    while (p < n - q && Math.abs(B.val[p*B.n + p + 1]) < TOL*(Math.abs(B.val[p*B.n + p]) + Math.abs(B.val[(p + 1)*B.n + (p + 1)]))) {
      p++;
    }

    if (q < n) {
      let DiagonalofB22isZero = -1;
      for (k = p; k < n - q; k++) {
        if (Math.abs(B.val[k*B.n + k]) < TOL) {
          DiagonalofB22isZero = k;
          break;
        }
      }
      if (DiagonalofB22isZero >= 0) {
        if (DiagonalofB22isZero < n - q - 1) {
          // Zero B(k,k+1) and entire row k...
          for (k = DiagonalofB22isZero + 1; k < n; k++) {

            cs = givens(B.val[k*B.n + k], B.val[DiagonalofB22isZero*B.n + k]);
            premulGivens(cs[0], cs[1], k, DiagonalofB22isZero, B);
            if (computeU)
              premulGivens(cs[0], cs[1], k, DiagonalofB22isZero, U);
          }
        } else {
          // Zero B(k-1,k) and entire column k...
          for (k = n - q - 2; k >= p; k--) {

            cs = givens(B.val[k*B.n + k], B.val[k*B.n + n - q - 1]);
            postmulGivens(cs[0], cs[1], k, n - q - 1, B);
            if (computeV)
              premulGivens(cs[0], cs[1], k, n - q - 1, Vt);
//							postmulGivens(cs[0],cs[1], j, n-q-1, V);
          }
        }
      } else {
        //B22 = get ( B, range(p , n - q ) , range (p , n-q ) );

        if (computeUV) {
          // UBV = GolubKahanSVDstep( B22, true ) ;
          // set ( U, range(p,n-q), [], mul(UBV.U, get(U, range(p,n-q), []) ) );
          // set ( Vt, range(p,n-q), [], mul(transpose(UBV.V), getRows(Vt, range(p,n-q)) ) );

          let GKstep = GolubKahanSVDstep(B, p, p, n - q - p, n - q - p, true);// this updates B22 inside B
          for (let kk = 0; kk < n - q - p - 1; kk++) {
            if (computeU)
              premulGivens(GKstep.csU[kk][0], GKstep.csU[kk][1], p + kk, p + kk + 1, U);
            if (computeV)
              premulGivens(GKstep.csV[kk][0], GKstep.csV[kk][1], p + kk, p + kk + 1, Vt); // premul because Vtransposed
          }
        } else {
          GolubKahanSVDstep(B, p, p, n - q - p, n - q - p);
        }
        //set ( B , range(p , n - q ) , range (p , n-q ), B22  );
      }
    }
    iter++;
  } while (q < n) ;

  if (computeUV) {

    if (computeV)
      V = transposeMatrix(Vt);

    // Correct sign of singular values:
    let s = diag(B);
    let signs = zeros(n);
    for (i = 0; i < n; i++) {
      if (s[i] < 0) {
        if (computeV)
          set(V, [], i, minus(get(V, [], i)));
        s[i] = -s[i];
      }
    }

    // Rearrange in decreasing order:
    let indexes = sort(s, true, true);
    if (computeV)
      V = get(V, [], indexes);
    if (computeU) {
      if (!thinU) {
        for (i = n; i < m; i++) {
          indexes.push(i);
        }
      }
      U = get(U, indexes, []);
    }

    let S;

    if (thinU)
      S = diag(s);
    else
      S = mat([diag(s), zeros(m - n, n)], true);

    let Ut = undefined;
    if (computeU)
      Ut = transpose(U);

    if (Atransposed) {
      if (thinU)
        return {"U": V, "S": S, "V": Ut, "s": s};
      else
        return {"U": V, "S": transpose(S), "V": Ut, "s": s};
    } else {
      return {"U": Ut, "S": S, "V": V, "s": s};
    }
  } else
    return sort(abs(diag(B)), true);
}

export function rank(A) {
  const s = svd(A);
  let rank = 0;
  let i;
  for (i = 0; i < s.length; i++) {
    if (s[i] > 1e-10)
      rank++;
  }

  return rank;
}

export function nullspace(A) {
  // Orthonormal basis for the null space of A
  const s = svd(A, "V");
  const n = A.n;

  let rank = 0;
  const TOL = 1e-8;
  while (rank < n && s.s[rank] > TOL) {
    rank++;
  }

  if (rank < n)
    return get(s.V, [], range(rank, n));
  else
    return zeros(n);

}

export function orth(A) {
  // Orthonormal basis for the range of A
  const s = svd(A, "thinU");
  const n = A.n;

  let rank = 0;
  const TOL = 1e-8;
  while (rank < n && s.s[rank] > TOL) {
    rank++;
  }

  return get(s.U, [], range(0, rank));

}

