import {type} from './lalolibbase.js';
import {
  zeros, ones, transpose, mul, add, issymmetric, isScalar, randn,
  randnScalar, chol, cholsolve, prodVector, prodMatrix, prodMatrixCols, prodMatrixRows,
  entrywisediv, entrywisemul, entrywisemulMatrix, entrywisemulVector, Matrix, addMatrices,
  appendRow, dot, diag, divMatrices, divMatrixScalar, divVectors, divVectorScalar, divScalarMatrix,
  divScalarVector
} from './linalg.js';

/**
 * Number of combinations of k items among n
 * @param{number}
 * @param{number}
 * @return{number}
 */
export function nchoosek(n, k) {
  if (k > n || k < 0 || n < 0)
    return 0;

  let i;
  let res = 1;
  for (i = n - k + 1; i <= n; i++) {
    res *= i;
  }
  for (i = 2; i <= k; i++) {
    res /= i;
  }
  return res;
}

//////////////////////////////////////////////
//  Multivariate Gaussian random vectors
//////////////////////////////////////////////
export function mvnrnd(mu, Sigma, N) {
  if (arguments.length < 3)
    let N = 1;

  let X = randn(N, mu.length);

  if (issymmetric(Sigma))
    let L = chol(Sigma);
  else
    let L = Sigma; // L directly provided instead of Sigma

  return add(mul(ones(N), transpose(mu)), mul(X, transpose(L)));
}

//////////////////////////////////////////////
// Generic class for Distributions
/////////////////////////////////////////////
export class Distribution {
  constructor(distrib, arg1, arg2) {

    if (arguments.length < 1) {
      error("Error in new Distribution(name): name is undefined.");
      return undefined;
    }

    if (typeof (distrib) === "string")
      distrib = eval(distrib);

    this.type = "Distribution:" + distrib.name;

    this.distribution = distrib.name;

    // Functions that depend on the distrib:
    this.construct = distrib.prototype.construct;

    this.estimate = distrib.prototype.estimate;
    this.sample = distrib.prototype.sample;
    this.pdf = distrib.prototype.pdf;
    if (distrib.prototype.pmf)
      this.pmf = distrib.prototype.pmf;

    if (distrib.prototype.logpdf)
      this.logpdf = distrib.prototype.logpdf;
    else
      this.logpdf = function (x) {
        return log(this.pdf(x));
      };

//	this.cdf = distrib.prototype.cdf; 

    // Initialization depending on distrib
    this.construct(arg1, arg2);
  }

  construct(params) {
    // Read params and create the required fields for a specific algorithm

  }

  pdf(x) {
    // return value of PDF at x
  }

  sample(N) {
    // Return N samples
  }

  estimate(X) {
    // Estimate dsitribution from the N-by-d matrix X
    // !! this function should also update this.mean and this.variance
  }


  info() {
    // Print information about the distribution

    let str = "{<br>";
    let i;
    let Functions = new Array();
    for (i in this) {
      switch (type(this[i])) {
        case "string":
        case "boolean":
        case "number":
          str += i + ": " + this[i] + "<br>";
          break;
        case "vector":
          str += i + ": " + printVector(this[i]) + "<br>";
          break;
        case "matrix":
          str += i + ": matrix of size " + this[i].m + "-by-" + this[i].n + "<br>";
          break;
        case "function":
          Functions.push(i);
          break;
        default:
          str += i + ": " + typeof (this[i]) + "<br>";
          break;
      }
    }
    str += "<i>Functions: " + Functions.join(", ") + "</i><br>";
    str += "}";
    return str;
  }
}


///////////////////////////////
///  Uniform 
///////////////////////////////
export class Uniform extends Distribution {
  constructor(params) {
    super(params);
  }


  construct(a, b) {
    // Read params and create the required fields for a Uniform distribution
    if (typeof (a) === "undefined") {
      // default to continuous uniform in [-1,1];
      this.isDiscrete = false;
      this.a = -1;
      this.b = 1;
      this.dimension = 1;

      this.px = 0.5;
      this.mean = 0;
      this.variance = 1/3;
      this.std = Math.sqrt(this.variance);
    } else {
      if (typeof (b) === "undefined") {
        this.isDiscrete = true;
        if (typeof (a) === "number")
          this.values = range(a);
        else
          this.values = a;
        this.dimension = 1;
        this.mean = (min(this.values) + max(this.values))/2;
        this.variance = (this.values.length*this.values.length - 1)/12;
        this.std = Math.sqrt(this.variance);
      } else {
        this.isDiscrete = false;
        this.a = a;
        this.b = b;
        this.dimension = size(a, 1);

        this.px = 1/prod(sub(b, a));
        this.mean = mul(0.5, add(a, b));
        let b_a = sub(b, a);
        this.variance = entrywisediv(entrywisemul(b_a, b_a), 12);
        this.std = sqrt(this.variance);
      }
    }
  }

  pdf(x) {
    // return value of PDF at x
    const tx = type(x);
    let p = undefined;
    if (this.isDiscrete) {
      function pdfscalar (s, values) {
        return (values.indexOf(s) < 0) ? 0 : (1/values.length);
      }


      if (tx === "number") {
        p = pdfscalar(x, this.values);
      } else if (tx === "vector") {
        p = zeros(x.length);
        for (let i = 0; i < x.length; i++) {
          p[i] = pdfscalar(x[i], this.values);
        }
      } else if (tx === "matrix") {
        p = zeros(x.m, x.n);
        for (let i = 0; i < x.m*x.n; i++) {
          p.val[i] = pdfscalar(x.val[i], this.values);
        }
      }
    } else {
      function pdfscalar(s, l, u, px) {
        return (s >= l && s <= u) ? px : 0;
      }


      if (tx === "number") {
        if (this.dimension === 1)
          p = pdfscalar(x, this.a, this.b, this.px);
      } else if (tx === "vector") {
        if (this.dimension === 1) {
          p = zeros(x.length);
          for (let i = 0; i < x.length; i++) {
            p[i] = pdfscalar(x[i], this.a, this.b, this.px);
          }
        } else if (this.dimension === x.length) {
          p = pdfscalar(x[0], this.a[0], this.b[0], this.px);
          let k = 1;
          while (k < x.length && p !== 0) {
            p *= pdfscalar(x[k], this.a[k], this.b[k], this.px);
            k++;
          }
        }
      } else if (tx === "matrix") {
        if (this.dimension === 1) {
          p = zeros(x.m, x.n);
          for (let i = 0; i < x.m*x.n; i++) {
            p.val[i] = pdfscalar(x.val[i], this.a, this.b, this.px);
          }
        } else if (this.dimension === x.n) {
          p = zeros(x.m);
          for (let i = 0; i < x.m; i++) {
            p[i] = pdfscalar(x.val[i*x.n], this.a[0], this.b[0], this.px);
            let k = 1;
            while (k < x.n && p[i] !== 0) {
              p[i] *= pdfscalar(x.val[i*x.n + k], this.a[k], this.b[k], this.px);
              k++;
            }
          }
        }
      }
    }
    return p;
  }

  sample(N) {
    // Return N samples
    if (typeof (N) === "undefined")
      let N = 1;

    if (this.isDiscrete) {
      let s = zeros(N);
      for (let i = 0; i < N; i++) {
        let r = Math.random();
        let k = 1;
        let n = this.values.length;
        while (r > k/n) {
          k++;
        }
        s[i] = this.values[k - 1];
      }
      if (N === 1)
        return s[0];
      else
        return s;
    } else {
      if (this.dimension === 1)
        return add(entrywisemul(this.b - this.a, rand(N)), this.a);
      else {
        return add(entrywisemul(outerprod(ones(N), sub(this.b, this.a)), rand(N, this.dimension)), outerprod(ones(N), this.a));
      }
    }
  }

  estimate(X) {
    // Estimate dsitribution from the N-by-d matrix X
    const tX = type(X);

    // Detect if X contains discrete or continuous values
    if (tX === "matrix")
      let x = X.val;
    else
      let x = X;

    let i = 0;
    while (i < x.length && Math.round(x[i]) === x[i]) {
      i++;
    }
    if (i < x.length)
      this.isDiscrete = false;
    else
      this.isDiscrete = true;

    // Estimate
    if (this.isDiscrete) {
      for (i = 0; i < x.length; i++) {
        let xi = Math.round(x[i]);
        if (this.values.indexOf(xi) < 0)
          this.values.push(xi);
      }
      this.dimension = 1;
      this.mean = (min(this.values) + max(this.values))/2;
      this.variance = (this.values.length*this.values.length - 1)/12;
      this.std = Math.sqrt(this.variance);
    } else {
      if (tX === "matrix") {
        this.a = min(X, 1).val;
        this.b = max(X).val;
        this.dimension = this.a.length;
      } else {
        this.a = minVector(X);
        this.b = maxVector(X);
        this.dimension = 1;
      }
      this.mean = mul(0.5, add(this.a, this.b));
      let b_a = sub(this.b, this.a);
      this.variance = entrywisediv(entrywisemul(b_a, b_a), 12);
      this.std = sqrt(this.variance);
      this.px = 1/prod(sub(this.b, this.a));
    }
    return this;
  }
}

///////////////////////////////
///  Gaussian 
/// (with independent components in multiple dimension)
///////////////////////////////
export class Guassian extends Distribution {
  constructor(params) {
    super(params);
  }


  construct(mean, variance) {
    // Read params and create the required fields for a specific algorithm
    if (typeof (mean) === "undefined")
      let mu = 1;

    else if (type(mean) === "matrix")
      let mu = mean.val;
    else
      let mu = mean;

    let dim = size(mu, 1);

    if (typeof (variance) === "undefined") {
      if (dim === 1)
        let variance = 1;
      else
        let variance = ones(dim);
    }

    this.mean = mu;
    this.variance = variance;
    this.std = sqrt(this.variance);
    this.dimension = dim;
  }

  pdf(x) {
    // return value of PDF at x
    if (this.dimension === 1) {
      if (typeof (x) === "number") {
        let diff = x - this.mean;
        return Math.exp(-diff*diff/(2*this.variance))/(this.std*Math.sqrt(2*Math.PI));
      } else {
        let diff = sub(x, this.mean);
        return entrywisediv(exp(entrywisediv(entrywisemul(diff, diff), -2*this.variance)), this.std*Math.sqrt(2*Math.PI));
      }
    } else {
      if (type(x) === "vector") {
        if (x.length !== this.dimension) {
          error("Error in Gaussian.pdf(x): x.length = " + x.length + " !== " + this.dimension + " = Gaussian.dimension.");
          return undefined;
        }
        let diff = subVectors(x, this.mean);
        let u = -0.5*dot(diff, divVectors(diff, this.variance));
        return Math.exp(u)/(Math.pow(2*Math.PI, 0.5*this.dimension)*Math.sqrt(prodVector(this.variance)));
      } else {
        if (x.n !== this.dimension) {
          error("Error in Gaussian.pdf(X): X.n = " + x.n + " !== " + this.dimension + " = Gaussian.dimension.");
          return undefined;
        }

        let p = zeros(x.m);
        let denominator = Math.pow(2*Math.PI, 0.5*this.dimension)*Math.sqrt(prodVector(this.variance));
        for (let i = 0; i < x.m; i++) {
          let diff = subVectors(x.row(i), this.mean);
          let u = -0.5*dot(diff, divVectors(diff, this.variance));
          p[i] = Math.exp(u)/denominator;
        }
        return p;
      }
    }
  }

  sample(N) {
    // Return N samples
    if (typeof (N) === "undefined")
      let N = 1;

    if (N === 1)
      let X = add(entrywisemul(this.std, randn(this.dimension)), this.mean);
    else {
      let N1 = ones(N);
      let X = add(entrywisemul(outerprod(N1, this.std), randn(N, this.dimension)), outerprod(N1, this.mean));
    }
    return X;
  }

  estimate(X) {
    // Estimate dsitribution from the N-by-d matrix X
    if (type(X) === "matrix") {
      this.mean = mean(X, 1).val;
      this.variance = variance(X, 1).val;
      this.std = undefined;
      this.dimension = X.n;
    } else {
      this.mean = mean(X);
      this.variance = variance(X);
      this.std = Math.sqrt(this.variance);
      this.dimension = 1;
    }
    return this;
  }
}

///////////////////////////////
///  Gaussian 
/// (with independent components in multiple dimension)
///////////////////////////////
export class mvGaussian extends Distribution {
  constructor(params) {
    super(params);
  }


  construct(mean, covariance) {
    // Read params and create the required fields for a specific algorithm
    if (typeof (mean) === "undefined")
      let mu = 1;

    else if (type(mean) === "matrix")
      let mu = mean.val;
    else
      let mu = mean;

    let dim = size(mu, 1);

    if (typeof (covariance) === "undefined") {
      if (dim === 1)
        let covariance = 1;
      else
        let covariance = eye(dim);
    }

    this.mean = mu;
    this.variance = covariance;
    this.dimension = dim;

    this.L = chol(this.variance);
    if (typeof (this.L) === "undefined")
      error("Error in new Distribution (mvGaussian, mu, Sigma): Sigma is not positive definite");

    this.det = det(this.variance);
  }

  pdf(x) {
    // return value of PDF at x
    if (this.dimension === 1) {
      if (typeof (x) === "number") {
        let diff = x - this.mean;
        return Math.exp(-diff*diff/(2*this.variance))/(Math.sqrt(2*this.variance*Math.PI));
      } else {
        let diff = sub(x, this.mean);
        return entrywisediv(exp(entrywisediv(entrywisemul(diff, diff), -2*this.variance)), Math.sqrt(2*this.variance*Math.PI));
      }
    } else {
      if (type(x) === "vector") {
        if (x.length !== this.dimension) {
          error("Error in mvGaussian.pdf(x): x.length = " + x.length + " !== " + this.dimension + " = mvGaussian.dimension.");
          return undefined;
        }
        let diff = subVectors(x, this.mean);
        let u = -0.5*dot(diff, cholsolve(this.L, diff));
        return Math.exp(u)/Math.sqrt(Math.pow(2*Math.PI, this.dimension)*this.det);
      } else {
        if (x.n !== this.dimension) {
          error("Error in Gaussian.pdf(X): X.n = " + x.n + " !== " + this.dimension + " = Gaussian.dimension.");
          return undefined;
        }

        let p = zeros(x.m);
        let denominator = Math.sqrt(Math.pow(2*Math.PI, this.dimension)*this.det);
        for (let i = 0; i < x.m; i++) {
          let diff = subVectors(x.row(i), this.mean);
          let u = -0.5*dot(diff, cholsolve(this.L, diff));
          p[i] = Math.exp(u)/denominator;
        }
        return p;
      }
    }
  }

  sample(N) {
    // Return N samples
    if (typeof (N) === "undefined")
      let N = 1;

    let X = add(mul(randn(N, this.dimension), transpose(this.L)), outerprod(ones(N), this.mean));

    if (N === 1)
      return X.val;
    else
      return X;
  }

  estimate(X) {
    // Estimate dsitribution from the N-by-d matrix X
    if (type(X) === "matrix") {
      this.mean = mean(X, 1).val;
      this.variance = cov(X);
      this.dimension = X.n;
      this.L = chol(this.variance);
      if (typeof (this.L) === "undefined")
        error("Error in mvGaussian.estimate(X): covariance estimate is not positive definite");

      this.det = det(this.variance);
      return this;
    } else {
      error("mvGaussian.estimate( X ) needs a matrix X");
    }
  }
}


///////////////////////////////
///  Bernoulli 
///////////////////////////////
export class Bernoulli extends Distribution {
  constructor(params) {
    super(params);
  }

  construct(mean) {
    // Read params and create the required fields for a specific algorithm
    if (typeof (mean) === "undefined")
      let mean = 0.5;

    let dim = size(mean, 1);

    this.mean = mean;
    this.variance = entrywisemul(mean, sub(1, mean));
    this.std = sqrt(this.variance);
    this.dimension = dim;
  }

  pmf(x) {
    return this.pdf(x);
  }

  pdf(x) {
    // return value of PDF at x
    const tx = type(x);

    function pdfscalar(s, mu) {
      if (s === 1)
        return mu;
      else if (s === 0)
        return (1 - mu);
      else
        return 0;
    }


    let p;

    if (this.dimension === 1) {
      if (tx === "number") {
        return pdfscalar(x, this.mean);
      } else if (tx === "vector") {
        p = zeros(x.length);
        for (let i = 0; i < x.length; i++) {
          p[i] = pdfscalar(x[i], this.mean);
        }
        return p;
      } else if (tx === "matrix") {
        let P = zeros(x.m, x.n);
        let mn = x.m*x.n;
        for (let k = 0; k < mn; k++) {
          P.val[k] = pdfscalar(x.val[k], this.mean);
        }
        return P;
      }
    } else {
      switch (tx) {
        case "vector":
          p = pdfscalar(x[0], this.mean[0]);
          for (let k = 1; k < this.dimension; k++) {
            p *= pdfscalar(x[k], this.mean[k]);
          }
          break;

        case "spvector":
          p = 1;
          for (let j = 0; j < x.ind[0]; j++) {
            p *= (1 - this.mean[j]);
          }
          for (let k = 0; k < x.val.length - 1; k++) {
            p *= this.mean[x.ind[k]];
            for (let j = x.ind[k] + 1; j < x.ind[k + 1]; j++) {
              p *= (1 - this.mean[j]);
            }
          }
          p *= this.mean[x.ind[k]];
          for (let j = x.ind[k] + 1; j < this.dimension; j++) {
            p *= (1 - this.mean[j]);
          }
          break;

        case "matrix":
          p = zeros(x.m);
          for (let i = 0; i < x.m; i++) {
            p[i] = pdfscalar(x.val[i*x.n], this.mean[0]);
            for (let k = 1; k < x.n; k++) {
              p[i] *= pdfscalar(x.val[i*x.n + k], this.mean[k]);
            }
          }
          break;
        case "spmatrix":
          p = ones(x.m);
          for (let i = 0; i < x.m; i++) {
            let xr = x.row(i);	// could be faster without this...
            for (let j = 0; j < xr.ind[0]; j++) {
              p[i] *= (1 - this.mean[j]);
            }
            for (let k = 0; k < xr.val.length - 1; k++) {
              p[i] *= this.mean[xr.ind[k]];
              for (let j = xr.ind[k] + 1; j < xr.ind[k + 1]; j++) {
                p[i] *= (1 - this.mean[j]);
              }
            }
            p[i] *= this.mean[xr.ind[k]];
            for (let j = xr.ind[k] + 1; j < this.dimension; j++) {
              p[i] *= (1 - this.mean[j]);
            }
          }
          break;
        default:
          p = undefined;
          break;
      }
      return p;
    }

  }

  logpdf() {
    return this.logpmf(...arguments);
  }

  logpmf(x) {
    // return value of logPDF at x
    const tx = type(x);
    let p;

    function logpdfscalar(s, mu) {
      if (s === 1)
        return Math.log(mu);
      else if (s === 0)
        return Math.log(1 - mu);
      else
        return -Infinity;
    }

    if (this.dimension === 1) {
      if (tx === "number") {
        return logpdfscalar(x, this.mean);
      } else if (tx === "vector") {
        p = zeros(x.length);
        for (let i = 0; i < x.length; i++) {
          p[i] = logpdfscalar(x[i], this.mean);
        }
        return p;
      } else if (tx === "matrix") {
        let P = zeros(x.m, x.n);
        let mn = x.m*x.n;
        for (let k = 0; k < mn; k++) {
          P.val[k] = logpdfscalar(x.val[k], this.mean);
        }
        return P;
      }
    } else {
      switch (tx) {
        case "vector": {
          p = 0;
          for (let k = 0; k < this.dimension; k++) {
            p += logpdfscalar(x[k], this.mean[k]);
          }
          break;
        }
        case "spvector": {
          p = 0;
          for (let j = 0; j < x.ind[0]; j++) {
            p += Math.log(1 - this.mean[j]);
          }
          for (let k = 0; k < x.val.length - 1; k++) {
            p += Math.log(this.mean[x.ind[k]]);
            for (let j = x.ind[k] + 1; j < x.ind[k + 1]; j++) {
              p += Math.log(1 - this.mean[j]);
            }
          }
          p += Math.log(this.mean[x.ind[k]]);
          for (let j = x.ind[k] + 1; j < this.dimension; j++) {
            p += Math.log(1 - this.mean[j]);
          }
          break;
        }

        case "matrix": {
          p = zeros(x.m);
          for (let i = 0; i < x.m; i++) {
            for (let k = 0; k < x.n; k++) {
              p[i] += logpdfscalar(x.val[i*x.n + k], this.mean[k]);
            }
          }
          break;
        }
        case "spmatrix": {
          p = zeros(x.m);
          for (let i = 0; i < x.m; i++) {
            let xr = x.row(i);	// could be faster without this...
            for (let j = 0; j < xr.ind[0]; j++) {
              p[i] += Math.log(1 - this.mean[j]);
            }
            for (let k = 0; k < xr.val.length - 1; k++) {
              p[i] += Math.log(this.mean[xr.ind[k]]);
              for (let j = xr.ind[k] + 1; j < xr.ind[k + 1]; j++) {
                p[i] += Math.log(1 - this.mean[j]);
              }
            }
            p[i] += Math.log(this.mean[xr.ind[k]]);
            for (let j = xr.ind[k] + 1; j < this.dimension; j++) {
              p[i] += Math.log(1 - this.mean[j]);
            }
          }
          break;
        }

        default:
          p = undefined;
          break;
      }
      return p;
    }

  }

  sample(N) {
    // Return N samples
    if (typeof (N) === "undefined" || N === 1) {
      return isLower(rand(this.dimension), this.mean);
    } else {
      return isLower(rand(N, this.dimension), outerprod(ones(N), this.mean));
    }
  }


  estimate(X) {
    // Estimate dsitribution from the N-by-d matrix X
    switch (type(X)) {
      case "matrix":
      case "spmatrix":
        this.mean = mean(X, 1).val;
        this.variance = entrywisemul(this.mean, sub(1, this.mean));
        this.std = sqrt(this.variance);
        this.dimension = X.n;
        break;
      case "vector":
      case "spvector":
        this.dimension = 1;
        this.mean = mean(X);
        this.variance = this.mean*(1 - this.mean);
        this.std = Math.sqrt(this.variance);
        break;
      default:
        error("Error in Bernoulli.estimate( X ): X must be a (sp)matrix or (sp)vector.");
        break;
    }
    return this;
  }
}

///////////////////////////////
///  Poisson 
///////////////////////////////
export class Poisson extends Distribution {
  constructor(params) {
    super(params);
  }

  construct(mean) {
    // Read params and create the required fields for a specific algorithm
    if (typeof (mean) === "undefined")
      let mean = 5;

    let dim = size(mean, 1);

    this.mean = mean;
    this.variance = this.mean;
    this.std = sqrt(this.variance);
    this.dimension = dim;
  }

  pmf(x) {
    return this.pdf(x);
  }

  pdf(x) {
    // return value of PDF at x
    const tx = type(x);

    function pdfscalar (s, lambda) {
      if (s < 0 || Math.round(s) !== s)
        return 0;
      else if (s === 0)
        return 1;
      else {
        let u = lambda;
        for (let k = 2; k <= s; k++) {
          u *= lambda/k;
        }
        return Math.exp(-lambda)*u;
      }
    }


    if (this.dimension === 1) {
      if (tx === "number") {
        return pdfscalar(x, this.mean);
      } else if (tx === "vector") {
        let p = zeros(x.length);
        for (let i = 0; i < x.length; i++) {
          p[i] = pdfscalar(x[i], this.mean);
        }
        return p;
      } else if (tx === "matrix") {
        let P = zeros(x.m, x.n);
        let mn = x.m*x.n;
        for (let k = 0; k < mn; k++) {
          P.val[k] = pdfscalar(x.val[k], this.mean);
        }
        return p;
      }
    } else {
      if (tx === "vector") {
        let p = pdfscalar(x[0], this.mean[0]);
        for (let k = 0; k < this.dimension; k++) {
          p *= pdfscalar(x[k], this.mean[k]);
        }

        return p;
      } else if (tx === "matrix") {
        let p = zeros(x.m);
        for (let i = 0; i < x.m; i++) {
          p[i] = pdfscalar(x.val[i*x.n], this.mean[0]);
          for (let k = 0; k < x.n; k++) {
            p[i] *= pdfscalar(x.val[i*x.n + k], this.mean[k]);
          }
        }
        return p;
      }
    }

  }

  sample(N) {
    // Return N samples
    let samplescalar
    (lambda)
    {
      let x = Math.random();
      let n = 0;
      const exp_lambda = Math.exp(-lambda);
      while (x > exp_lambda) {
        x *= Math.random();
        n++;
      }
      return n;
    }
    ;

    if (typeof (N) === "undefined" || N === 1) {
      if (this.dimension === 1)
        return samplescalar(this.mean);
      else {
        let s = zeros(this.dimension);
        for (let k = 0; k < this.dimension; k++) {
          s[k] = samplescalar(this.mean[k]);
        }
        return s;
      }
    } else {
      if (this.dimension === 1) {
        let S = zeros(N);
        for (let i = 0; i < N; i++) {
          S[i] = samplescalar(this.mean);
        }
        return S;
      } else {
        let S = zeros(N, this.dimension);
        for (let i = 0; i < N; i++) {
          for (let k = 0; k < this.dimension; k++) {
            S[i*this.dimension + k] = samplescalar(this.mean[k]);
          }
        }
        return S;
      }
    }
  }


  estimate(X) {
    // Estimate dsitribution from the N-by-d matrix X
    if (type(X) === "matrix") {
      this.mean = mean(X, 1).val;
      this.variance = this.mean;
      this.std = sqrt(this.variance);
      this.dimension = X.n;
    } else { // X is a vector samples
      this.dimension = 1;
      this.mean = mean(X);
      this.variance = this.mean;
      this.std = Math.sqrt(this.variance);
    }
    return this;
  }
}
