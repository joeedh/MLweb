const Complex_I = new Complex(0, 1);

import {Matrix} from './linalg.js';

/**
 * @constructor
 * @struct
 */
export class Complex {
  constructor(a, b, polar) {
    /** @const */ this.type = "Complex";

    if (typeof (a) === "undefined") {
      this.re = 0.0;
      this.im = 0.0;
    } else if (a instanceof Complex) {
      this.re = a.re;
      this.im = a.im;
    } else if (typeof (a) === "number" && !polar) {
      this.re = a;
      this.im = b;
    } else {
      this.re = a*Math.cos(b);
      this.im = a*Math.sin(b);
    }
  }

  toString() {
    return this.re + (this.im >= 0 ? " + " : " - ") + Math.abs(this.im) + "i";
  }

  info() {
    return this.re + (this.im >= 0 ? " + " : " - ") + Math.abs(this.im) + "i";
  }
}

/**
 * @param {Complex}
 * @param {Complex}
 * @return {Complex}
 */
export function addComplex(a, b) {
  let z = new Complex(a);
  z.re += b.re;
  z.im += b.im;
  return z;
}

/**
 * @param {Complex}
 * @param {number}
 * @return {Complex}
 */
export function addComplexReal(a, b) {
  let z = new Complex(a);
  z.re += b;
  return z;
}

/**
 * @param {Complex}
 * @param {Complex}
 * @return {Complex}
 */
export function subComplex(a, b) {
  let z = new Complex(a);
  z.re -= b.re;
  z.im -= b.im;
  return z;
}

/**
 * @param {Complex}
 * @return {Complex}
 */
export function minusComplex(a) {
  return new Complex(-a.re, -a.im);
}

export function mulComplex(a, b) {
  return new Complex(a.re*b.re - a.im*b.im, a.im*b.re + a.re*b.im);
}

export function mulComplexReal(a, b) {
  return new Complex(a.re*b, a.im*b);
}

export function divComplex(a, b) {
  let denom = b.re*b.re + b.im*b.im;
  return new Complex((a.re*b.re + a.im*b.im)/denom, (a.im*b.re - a.re*b.im)/denom);
}

export function conj(z) {
  if (z instanceof Complex)
    return new Complex(z.re, -z.im);
  else if (z instanceof ComplexVector) {
    let r = new ComplexVector(z);
    for (let i = 0; i < z.length; i++) {
      r.im[i] = -r.im[i];
    }
    return r;
  } else if (z instanceof ComplexMatrix) {
    let r = new ComplexMatrix(z);
    for (let i = 0; i < z.length; i++) {
      r.im[i] = -r.im[i];
    }
    return r;
  } else
    return new Complex(z);	// for a real
}

export function modulus(z) {
  if (z instanceof Complex)
    return Math.sqrt(z.re*z.re + z.im*z.im);
  else if (z instanceof ComplexVector)
    return sqrt(addVectors(entrywisemulVector(z.re, z.re), entrywisemulVector(z.im, z.im)));
  else if (z instanceof ComplexVector)
    return new Matrix(z.m, z.n, sqrt(addVectors(entrywisemulVector(z.re, z.re), entrywisemulVector(z.im, z.im)), true));
}

let absComplex = modulus;

export function expComplex(z) {
  return new Complex(Math.exp(z.re), z.im, true);
}


/**
 * @constructor
 * @struct
 */
export class ComplexVector {
  constructor(a, b, dontcopy) {
    /** @const */ this.type = "ComplexVector";

    if (arguments.length === 0) {
      // dummy call, probably in renewObject
      // while loading data from a file
    } else if (a instanceof ComplexVector) {
      /** @const */ this.length = a.length;
      this.re = vectorCopy(a.re);
      this.im = vectorCopy(a.im);
    } else if (typeof (a) === "number") {
      /** @const */ this.length = a;
      this.re = new Float64Array(a);
      this.im = new Float64Array(a);
    } else if (a instanceof Float64Array && b instanceof Float64Array) {
      /** @const */ this.length = a.length;
      if (typeof (dontcopy) === "undefined" || !dontcopy) {
        this.re = vectorCopy(a);
        this.im = vectorCopy(b);
      } else {
        this.re = a;
        this.im = b;
      }
    } else {
      error("Bad arguments to new ComplexVector()");
    }
  }

  toString() {
    return "[" + this.type + " of size " + this.length + "]";
  }

  get(i) {
    return new Complex(this.re[i], this.im[i]);
  }

  set(i, z) {
    if (typeof (z) === "number") {
      this.re[i] = z;
      this.im[i] = 0;
    } else {
      this.re[i] = z.re;
      this.im[i] = z.im;
    }
  }

  getSubVector(rowsrange) {
    const n = rowsrange.length;
    let res = new ComplexVector(n);
    for (let i = 0; i < n; i++) {
      res.re[i] = this.re[rowsrange[i]];
      res.im[i] = this.im[rowsrange[i]];
    }
    return res;
  }

  setVectorScalar(rowsrange, B) {
    let i;
    for (i = 0; i < rowsrange.length; i++) {
      this.set(rowsrange[i], B);
    }
  }

  setVectorVector(rowsrange, B) {
    let i;
    for (i = 0; i < rowsrange.length; i++) {
      this.set(rowsrange[i], B[i]);
    }
  }
}

/**
 * @constructor
 * @struct
 */
export class ComplexMatrix {
  constructor(a, b, values, valuesimag) {
    /** @const */ this.type = "ComplexMatrix";

    if (arguments.length === 0) {
      // dummy call, probably in renewObject
      // while loading data from a file
    } else if (a instanceof ComplexMatrix) {
      /** @const */ this.length = a.length;
      /** @const */ this.m = a.m;
      /** @const */ this.n = a.n;
      /** @const */ this.size = [a.m, a.n];
      this.re = vectorCopy(a.re);
      this.im = vectorCopy(a.im);
    } else if (typeof (a) === "number" && typeof (b) === "number") {
      /** @const */ this.length = a;
      /** @const */ this.m = a;
      /** @const */ this.n = b;
      /** @const */ this.size = [a, b];
      if (typeof (values) === "undefined") {
        this.re = new Float64Array(a*b);
        this.im = new Float64Array(a*b);
      } else if (values instanceof ComplexVector) {
        this.re = vectorCopy(values.re);
        this.im = vectorCopy(values.im);
      } else if (values instanceof Float64Array && typeof (valuesimag) !== "undefined" && valuesimag instanceof Float64Array) {
        this.re = values;
        this.im = valuesimag;	// !! no copy!
      }
    } else if (a instanceof Matrix && b instanceof Matrix) {
      /** @const */ this.length = a.length;
      /** @const */ this.m = a.m;
      /** @const */ this.n = a.n;
      /** @const */ this.size = [a.m, a.n];
      this.re = vectorCopy(a.val);
      this.im = vectorCopy(b.val);
    } else
      error("Bad arguments to new ComplexMatrix()");
  }


  toString() {
    return "[" + this.type + " of size " + this.m + " x " + this.n + "]";
  }

  get(i, j) {
    return new Complex(this.re[i*this.n + j], this.im[i*this.n + j]);
  }

  set(i, j, z) {
    if (typeof (z) === "number") {
      this.re[i*this.n + j] = z;
      this.im[i*this.n + j] = 0;
    } else {
      this.re[i*this.n + j] = z.re;
      this.im[i*this.n + j] = z.im;
    }
  }

  /**
   * @param {MatrixComplex}
   */
  transpose() {
    // simple Transpose without conjugate
    const m = A.m;
    const n = A.n;
    if (m > 1) {
      let i;
      let j;
      let res = new ComplexMatrix(n, m);
      let Aj = 0;
      for (j = 0; j < m; j++) {
        let ri = 0;
        for (i = 0; i < n; i++) {
          res.re[ri + j] = A.re[Aj + i];
          res.im[ri + j] = A.im[Aj + i];
          ri += m;
        }
        Aj += n;
      }
      return res;
    } else {
      return new ComplexVector(A.re, A.im);
    }
  }
}


export function real(z) {
  if (z instanceof Complex)
    return z.re;
  else if (z instanceof ComplexVector)
    return vectorCopy(z.re);
  else if (z instanceof ComplexMatrix)
    return new Matrix(z.m, z.n, z.re);
  else
    return copy(z);
}

export function imag(z) {
  if (z instanceof Complex)
    return z.im;
  else if (z instanceof ComplexVector)
    return vectorCopy(z.im);
  else if (z instanceof ComplexMatrix)
    return new Matrix(z.m, z.n, z.im);
  else
    return 0;
}

/**
 * @param {MatrixComplex}
 */
export function transposeComplexMatrix(A) {
  // Hermitian transpose = conjugate transpose
  const m = A.m;
  const n = A.n;
  if (m > 1) {
    let i;
    let j;
    let res = new ComplexMatrix(n, m);
    let Aj = 0;
    for (j = 0; j < m; j++) {
      let ri = 0;
      for (i = 0; i < n; i++) {
        res.re[ri + j] = A.re[Aj + i];
        res.im[ri + j] = -A.im[Aj + i];
        ri += m;
      }
      Aj += n;
    }
    return res;
  } else {
    return new ComplexVector(A.re, minusVector(A.im));
  }
}


/**
 * @param {ComplexVector}
 * @param {ComplexVector}
 * @return {ComplexVector}
 */
export function addComplexVectors(a, b) {
  let z = new ComplexVector(a);
  const n = a.length;
  for (let i = 0; i < n; i++) {
    z.re[i] += b.re[i];
    z.im[i] += b.im[i];
  }
  return z;
}

/**
 * @param {ComplexVector}
 * @param {ComplexVector}
 * @return {ComplexVector}
 */
export function subComplexVectors(a, b) {
  let z = new ComplexVector(a);
  const n = a.length;
  for (let i = 0; i < n; i++) {
    z.re[i] -= b.re[i];
    z.im[i] -= b.im[i];
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @param {ComplexMatrix}
 * @return {ComplexMatrix}
 */
export function addComplexMatrices(a, b) {
  let z = new ComplexMatrix(a);
  const mn = a.m*a.n;
  for (let i = 0; i < mn; i++) {
    z.re[i] += b.re[i];
    z.im[i] += b.im[i];
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @param {ComplexMatrix}
 * @return {ComplexMatrix}
 */
export function subComplexMatrices(a, b) {
  let z = new ComplexMatrix(a);
  const mn = a.m*a.n;
  for (let i = 0; i < mn; i++) {
    z.re[i] -= b.re[i];
    z.im[i] -= b.im[i];
  }
  return z;
}

/**
 * @param {ComplexVector}
 * @param {Float64Array}
 * @return {ComplexVector}
 */
export function addComplexVectorVector(a, b) {
  let z = new ComplexVector(a);
  const n = a.length;
  for (let i = 0; i < n; i++) {
    z.re[i] += b[i];
  }
  return z;
}

/**
 * @param {ComplexVector}
 * @param {Float64Array}
 * @return {ComplexVector}
 */
export function subComplexVectorVector(a, b) {
  let z = new ComplexVector(a);
  const n = a.length;
  for (let i = 0; i < n; i++) {
    z.re[i] -= b[i];
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @param {Matrix}
 * @return {ComplexMatrix}
 */
export function addComplexMatrixMatrix(a, b) {
  let z = new ComplexMatrix(a);
  const n = a.m*a.n;
  for (let i = 0; i < n; i++) {
    z.re[i] += b.val[i];
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @param {Matrix}
 * @return {ComplexMatrix}
 */
export function subComplexMatrixMatrix(a, b) {
  let z = new ComplexMatrix(a);
  const n = a.m*a.n;
  for (let i = 0; i < n; i++) {
    z.re[i] -= b.val[i];
  }
  return z;
}

/**
 * @param {number}
 * @param {ComplexVector}
 * @return {ComplexVector}
 */
export function addScalarComplexVector(a, b) {
  let z = new ComplexVector(b);
  const n = b.length;
  for (let i = 0; i < n; i++) {
    z.re[i] += a;
  }
  return z;
}

/**
 * @param {number}
 * @param {ComplexVector}
 * @return {ComplexVector}
 */
export function subScalarComplexVector(a, b) {
  let z = minusComplexVector(b);
  const n = b.length;
  for (let i = 0; i < n; i++) {
    z.re[i] += a;
  }
  return z;
}

/**
 * @param {number}
 * @param {ComplexMatrix}
 * @return {ComplexMatrix}
 */
export function addScalarComplexMatrix(a, b) {
  let z = new ComplexMatrix(b);
  const n = b.m*b.n;
  for (let i = 0; i < n; i++) {
    z.re[i] += a;
  }
  return z;
}


/**
 * @param {ComplexVector}
 * @param {ComplexVector}
 * @return {ComplexVector}
 */
export function entrywisemulComplexVectors(a, b) {
  const n = a.length;
  let z = new ComplexVector(n);
  for (let i = 0; i < n; i++) {
    z.re[i] = a.re[i]*b.re[i] - a.im[i]*b.im[i];
    z.im[i] = a.im[i]*b.re[i] + a.re[i]*b.im[i];
  }
  return z;
}

/**
 * @param {ComplexVector}
 * @param {ComplexVector}
 * @return {ComplexVector}
 */
export function entrywisedivComplexVectors(a, b) {
  const n = a.length;
  let z = new ComplexVector(n);
  for (let i = 0; i < n; i++) {
    let bre = b.re[i];
    let bim = b.im[i];
    let denom = bre*bre + bim*bim;
    z.re[i] = (a.re[i]*bre + a.im[i]*bim)/denom;
    z.im[i] = (a.im[i]*bre - a.re[i]*bim)/denom;
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @param {ComplexMatrix}
 * @return {ComplexMatrix}
 */
export function entrywisemulComplexMatrices(a, b) {
  const n = a.m*a.n;
  let z = new ComplexMatrix(a.m, a.n);
  for (let i = 0; i < n; i++) {
    z.re[i] = a.re[i]*b.re[i] - a.im[i]*b.im[i];
    z.im[i] = a.im[i]*b.re[i] + a.re[i]*b.im[i];
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @param {ComplexMatrix}
 * @return {ComplexMatrix}
 */
export function entrywisedivComplexMatrices(a, b) {
  const n = a.m*a.n;
  let z = new ComplexMatrix(a.m, a.n);
  for (let i = 0; i < n; i++) {
    let bre = b.re[i];
    let bim = b.im[i];
    let denom = bre*bre + bim*bim;
    z.re[i] = (a.re[i]*bre + a.im[i]*bim)/denom;
    z.im[i] = (a.im[i]*bre - a.re[i]*bim)/denom;
  }
  return z;
}

/**
 * @param {ComplexVector}
 * @param {Float64Array}
 * @return {ComplexVector}
 */
export function entrywisemulComplexVectorVector(a, b) {
  const n = a.length;
  let z = new ComplexVector(n);
  for (let i = 0; i < n; i++) {
    z.re[i] = a.re[i]*b[i];
    z.im[i] = a.im[i]*b[i];
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @param {Matrix}
 * @return {ComplexMatrix}
 */
export function entrywisemulComplexMatrixMatrix(a, b) {
  const n = a.m*a.n;
  let z = new ComplexMatrix(a.m, a.n);
  for (let i = 0; i < n; i++) {
    z.re[i] = a.re[i]*b.val[i];
    z.im[i] = a.im[i]*b.val[i];
  }
  return z;
}

/**
 * @param {ComplexVector}
 * @return {ComplexVector}
 */
export function minusComplexVector(a) {
  const n = a.length;
  let z = new ComplexVector(n);
  for (let i = 0; i < n; i++) {
    z.re[i] = -a.re[i];
    z.im[i] = -a.im[i];
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @return {ComplexMatrix}
 */
export function minusComplexMatrix(a) {
  let z = new ComplexMatrix(a.m, a.n);
  const n = a.m*a.n;
  for (let i = 0; i < n; i++) {
    z.re[i] = -a.re[i];
    z.im[i] = -a.im[i];
  }
  return z;
}

/**
 * @param {ComplexVector}
 * @return {number}
 */
export function sumComplexVector(a) {
  let z = new Complex();
  const n = a.length;
  for (let i = 0; i < n; i++) {
    z.re += a.re[i];
    z.im += a.im[i];
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @return {number}
 */
export function sumComplexMatrix(a) {
  let z = new Complex();
  const n = a.m*a.n;
  for (let i = 0; i < n; i++) {
    z.re += a.re[i];
    z.im += a.im[i];
  }
  return z;
}

/**
 * @param {ComplexVector}
 * @return {number}
 */
export function norm1ComplexVector(a) {
  let r = 0.0;
  const n = a.length;
  for (let i = 0; i < n; i++) {
    r += Math.sqrt(a.re[i]*a.re[i] + a.im[i]*a.im[i]);
  }
  return r;
}

/**
 * @param {ComplexVector}
 * @return {number}
 */
export function norm2ComplexVector(a) {
  let r = 0.0;
  const n = a.length;
  for (let i = 0; i < n; i++) {
    r += a.re[i]*a.re[i] + a.im[i]*a.im[i];
  }
  return Math.sqrt(r);
}

/**
 * @param {ComplexMatrix}
 * @return {number}
 */
export function normFroComplexMatrix(a) {
  let r = 0.0;
  const n = a.m*a.n;
  for (let i = 0; i < n; i++) {
    r += a.re[i]*a.re[i] + a.im[i]*a.im[i];
  }
  return Math.sqrt(r);
}

/**
 * @param {ComplexVector}
 * @param {ComplexVector}
 * @return {Complex}
 */
export function dotComplexVectors(a, b) {
  // = b^H a = conj(b)^T a
  let z = new Complex();
  const n = a.length;
  for (let i = 0; i < n; i++) {
    z.re += a.re[i]*b.re[i] + a.im[i]*b.im[i];
    z.im += a.im[i]*b.re[i] - a.re[i]*b.im[i]
  }
  return z;
}

/**
 * @param {ComplexVector}
 * @param {Float64Array}
 * @return {Complex}
 */
export function dotComplexVectorVector(a, b) {
  // = b^T a
  let z = new Complex();
  const n = a.length;
  for (let i = 0; i < n; i++) {
    z.re += a.re[i]*b[i];
    z.im += a.im[i]*b[i];
  }
  return z;
}

/**
 * @param {number}
 * @param {ComplexVector}
 * @return {ComplexVector}
 */
export function mulScalarComplexVector(a, b) {
  let re = mulScalarVector(a, b.re);
  let im = mulScalarVector(a, b.im);
  return new ComplexVector(re, im, true);
}

/**
 * @param {Complex}
 * @param {ComplexVector}
 * @return {ComplexVector}
 */
export function mulComplexComplexVector(a, b) {
  const n = b.length;
  let z = new ComplexVector(n);
  let are = a.re;
  let aim = a.im;
  for (let i = 0; i < n; i++) {
    z.re[i] = are*b.re[i] - aim*b.im[i];
    z.im[i] = aim*b.re[i] + are*b.im[i];
  }
  return z;
}

/**
 * @param {Complex}
 * @param {Float64Array}
 * @return {ComplexVector}
 */
export function mulComplexVector(a, b) {
  const n = b.length;
  let z = new ComplexVector(n);
  let are = a.re;
  let aim = a.im;
  for (let i = 0; i < n; i++) {
    z.re[i] = are*b[i];
    z.im[i] = aim*b[i];
  }
  return z;
}

/**
 * @param {number}
 * @param {ComplexMatrix}
 * @return {ComplexMatrix}
 */
export function mulScalarComplexMatrix(a, b) {
  let re = mulScalarVector(a, b.re);
  let im = mulScalarVector(a, b.im);
  return new ComplexMatrix(b.m, b.n, re, im);
}

/**
 * @param {Complex}
 * @param {ComplexMatrix}
 * @return {ComplexMatrix}
 */
export function mulComplexComplexMatrix(a, b) {
  const n = b.m*b.n;
  let z = new ComplexMatrix(b.m, b.n);
  let are = a.re;
  let aim = a.im;
  for (let i = 0; i < n; i++) {
    z.re[i] = are*b.re[i] - aim*b.im[i];
    z.im[i] = aim*b.re[i] + are*b.im[i];
  }
  return z;
}

/**
 * @param {Complex}
 * @param {Matrix}
 * @return {ComplexMatrix}
 */
export function mulComplexMatrix(a, b) {
  const n = b.m*b.n;
  let z = new ComplexMatrix(b.m, b.n);
  let are = a.re;
  let aim = a.im;
  for (let i = 0; i < n; i++) {
    z.re[i] = are*b.val[i];
    z.im[i] = aim*b.val[i];
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @param {Float64Array}
 * @return {ComplexVector}
 */
export function mulComplexMatrixVector(a, b) {
  const m = a.m;
  const n = a.n;
  let z = new ComplexVector(m);
  let ai = 0;
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      z.re[i] += a.re[ai + j]*b[j];
      z.im[i] += a.im[ai + j]*b[j];
    }
    ai += n;
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @param {ComplexVector}
 * @return {ComplexVector}
 */
export function mulComplexMatrixComplexVector(a, b) {
  const m = a.m;
  const n = a.n;
  let z = new ComplexVector(m);
  let ai = 0;
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < n; j++) {
      z.re[i] += a.re[ai + j]*b.re[j] - a.im[ai + j]*b.im[j];
      z.im[i] += a.im[ai + j]*b.re[j] + a.re[ai + j]*b.im[j];
    }
    ai += n;
  }
  return z;
}

/**
 * @param {ComplexMatrix}
 * @param {ComplexMatrix}
 * @return {ComplexMatrix}
 */
export function mulComplexMatrices(A, B) {
  const m = A.length;
  const n = B.n;
  const n2 = B.length;

  let Are = A.re;
  let Aim = A.im;
  let Bre = B.re;
  let Bim = B.im;

  let Cre = new Float64Array(m*n);
  let Cim = new Float64Array(m*n);
  let aik;
  let Aik = 0;
  let Ci = 0;
  for (let i = 0; i < m; i++) {
    let bj = 0;
    for (let k = 0; k < n2; k++) {
      let aikre = Are[Aik];
      let aikim = Aim[Aik];

      for (let j = 0; j < n; j++) {
        Cre[Ci + j] += aikre*Bre[bj] - aikim*Bim[bj];
        Cim[Ci + j] += aikre*Bim[bj] + aikim*Bre[bj];
        bj++;
      }
      Aik++;
    }
    Ci += n;
  }
  return new ComplexMatrix(m, n, Cre, Cim);
}

/**
 * @param {ComplexMatrix}
 * @param {Matrix}
 * @return {ComplexMatrix}
 */
export function mulComplexMatrixMatrix(A, B) {
  const m = A.m;
  const n = B.n;
  const n2 = B.m;

  let Are = A.re;
  let Aim = A.im;
  let Bre = B.val;

  let Cre = new Float64Array(m*n);
  let Cim = new Float64Array(m*n);
  let aik;
  let Aik = 0;
  let Ci = 0;
  for (let i = 0; i < m; i++) {
    let bj = 0;
    for (let k = 0; k < n2; k++) {
      let aikre = Are[Aik];
      let aikim = Aim[Aik];

      for (let j = 0; j < n; j++) {
        Cre[Ci + j] += aikre*Bre[bj];
        Cim[Ci + j] += aikim*Bre[bj];
        bj++;
      }
      Aik++;
    }
    Ci += n;
  }
  return new ComplexMatrix(m, n, Cre, Cim);
}


/**
 * @param {Float64Array|ComplexVector}
 * @return {ComplexVector}
 */
export function fft(x) {
  const n = x.length;
  const s = Math.log2(n);
  const m = n/2;

  if (s%1 !== 0) {
    error("fft(x) only implemented for x.length = 2^m. Use dft(x) instead.");
    return undefined;
  }

  let X = new ComplexVector(x, zeros(n));

  // bit reversal:
  let j = 0;
  for (let i = 0; i < n - 1; i++) {
    if (i < j) {
      // swap(X[i], X[j])
      let Xi = X.re[i];
      X.re[i] = X.re[j];
      X.re[j] = Xi;
      Xi = X.im[i];
      X.im[i] = X.im[j];
      X.im[j] = Xi;
    }

    let k = m;
    while (k <= j) {
      j -= k;
      k /= 2;
    }
    j += k;
  }

  // FFT:
  let l2 = 1;
  let c = new Complex(-1, 0);
  let u = new Complex();
  for (let l = 0; l < s; l++) {
    let l1 = l2;
    l2 *= 2;
    u.re = 1;
    u.im = 0;
    for (let j = 0; j < l1; j++) {
      for (let i = j; i < n; i += l2) {
        let i1 = i + l1;
        //let t1 = mulComplex(u, X.get(i1) );
        let t1re = u.re*X.re[i1] - u.im*X.im[i1]; // t1 = u * X[i1]
        let t1im = u.im*X.re[i1] + u.re*X.im[i1];

        X.re[i1] = X.re[i] - t1re;
        X.im[i1] = X.im[i] - t1im;

        X.re[i] += t1re;
        X.im[i] += t1im;
      }

      u = mulComplex(u, c);
    }

    c.im = -Math.sqrt((1.0 - c.re)/2.0);
    c.re = Math.sqrt((1.0 + c.re)/2.0);
  }
  return X;
}

/**
 * @param {ComplexVector}
 * @return {ComplexVector|Float64Array}
 */
export function ifft(x) {
  const n = x.length;
  const s = Math.log2(n);
  const m = n/2;

  if (s%1 !== 0) {
    error("ifft(x) only implemented for x.length = 2^m. Use idft(x) instead.");
    return undefined;
  }


  let X = new ComplexVector(x, zeros(n));

  // bit reversal:
  let j = 0;
  for (let i = 0; i < n - 1; i++) {
    if (i < j) {
      // swap(X[i], X[j])
      let Xi = X.re[i];
      X.re[i] = X.re[j];
      X.re[j] = Xi;
      Xi = X.im[i];
      X.im[i] = X.im[j];
      X.im[j] = Xi;
    }

    let k = m;
    while (k <= j) {
      j -= k;
      k /= 2;
    }
    j += k;
  }

  // iFFT:
  let l2 = 1;
  let c = new Complex(-1, 0);
  let u = new Complex();
  for (let l = 0; l < s; l++) {
    let l1 = l2;
    l2 *= 2;
    u.re = 1;
    u.im = 0;
    for (let j = 0; j < l1; j++) {
      for (let i = j; i < n; i += l2) {
        let i1 = i + l1;
        //let t1 = mulComplex(u, X.get(i1) );
        let t1re = u.re*X.re[i1] - u.im*X.im[i1]; // t1 = u * X[i1]
        let t1im = u.im*X.re[i1] + u.re*X.im[i1];

        X.re[i1] = X.re[i] - t1re;
        X.im[i1] = X.im[i] - t1im;

        X.re[i] += t1re;
        X.im[i] += t1im;
      }

      u = mulComplex(u, c);
    }

    c.im = Math.sqrt((1.0 - c.re)/2.0);
    c.re = Math.sqrt((1.0 + c.re)/2.0);
  }

  let isComplex = false;
  for (let i = 0; i < n; i++) {
    X.re[i] /= n;
    X.im[i] /= n;
    if (Math.abs(X.im[i]) > 1e-6)
      isComplex = true;
  }
  if (isComplex)
    return X;
  else
    return X.re;
}

export function dft(x) {
  // DFT of a real signal
  if (typeof (x) === "number")
    return new Complex(x, 0);

  const n = x.length;
  if (n === 1)
    return new Complex(x[0], 0);
  else if (Math.log2(n)%1 === 0)
    return fft(x);
  else {
    let X = new ComplexVector(n);
    let thet = 0.0;
    for (let i = 0; i < n; i++) {
      let theta = 0.0;
      for (let t = 0; t < n; t++) {
        // theta = -2 pi i * t / n;
        X.re[i] += x[t]*Math.cos(theta);
        X.im[i] += x[t]*Math.sin(theta);
        theta += thet;
      }
      thet -= 2*Math.PI/n;
    }
    return X;
  }
}

export function idft(X) {
  // Only recovers real part
  /*
importScripts("src/experimental/complex.js")
t = 0:512
x = sin(t)
X = dft(x)
plot(modulus(X))
s = idft(X)
plot(s)	
  */
  if (!(X instanceof ComplexVector)) {
    if (X instanceof Complex)
      return X.re;
    else if (typeof (X) === "number")
      return X;
    else if (X instanceof Float64Array)
      return idft(new ComplexVector(X, zeros(X.length), true));
    else
      return undefined;
  }
  const n = X.length;
  if (n === 1)
    return X.re[0];
  else if (Math.log2(n)%1 === 0)
    return ifft(X);
  else {
    let x = new Float64Array(n);
    let thet = 0.0;
    for (let t = 0; t < n; t++) {
      let theta = 0.0;
      let re = 0.0;
      //let im = 0.0;
      for (let i = 0; i < n; i++) {
        // theta = 2 pi i * t / n;
        re += X.re[i]*Math.cos(theta) - X.im[i]*Math.sin(theta);
        // im += X[i].im * Math.sin(theta) + X[i].re * Math.cos(theta); // not used for real signals
        theta += thet;
      }
      x[t] = re/n;
      thet += 2*Math.PI/n;
    }
    return x;
  }
}

export function spectrum(x) {
  if (x instanceof Float64Array) {
    return absComplex(dft(x));
  } else
    return undefined;
}
