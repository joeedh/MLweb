
//////////////////////////
//// Cross-browser compatibility
///////////////////////////

if( typeof(console) == "undefined" ) {
  // for Safari
  var console = {log: function ( ) { } };
}

if( typeof(Math.sign) == "undefined" ) {
  // for IE, Safari
  Math.sign = function ( x ) { return ( x>=0 ? (x==0 ? 0 : 1) : -1 ) ;}
}

if (Array.prototype.remove === undefined) {
  Array.prototype.remove = function(item, no_error) {
    let i = this.indexOf(item);

    if (i < 0) {
      if (no_error) {
        console.warn("Missing item", item);
      } else {
        throw new Error("Missing item");
      }
    }

    while (i < this.length-1) {
      this[i] = this[i+1];
      i++;
    }

    this[i] = undefined;
    this.length--;

    return this;
  }
}

Math.fract = f => f - Math.floor(f);
Math.tent = f => 1.0 - Math.abs(Math.fract(f)-0.5)*2.0;

if (Set.prototype.map === undefined) {
  Set.prototype.map = function(cb) {
    let ret = new Set();

    for (let item of this) {
      ret.add(cb(item));
    }

    return ret;
  }
}

if (Set.prototype.filter === undefined) {
  Set.prototype.filter = function(cb) {
    let ret = new Set();

    for (let item of this) {
      if (cb(item)) {
        ret.add(item);
      }
    }

    return ret;
  }
}

export function makeFastArrayIters() {
  let arrayiters = new Array(1024);
  arrayiters.cur = 0;

  class ArrayIter {
    constructor() {
      this.array = undefined;
      this.i = 0;
      this.done = false;
      this.ret = {value: undefined, done: true};
    }

    reset(array) {
      this.array = array;
      this.i = 0;
      this.done = false;
      this.ret.done = false;

      return this;
    }

    [Symbol.iterator]() {
      return this;
    }

    next() {
      if (this.i >= this.array.length) {
        return this.finish();
      }

      this.ret.value = this.array[this.i++];
      return this.ret;
    }

    finish() {
      if (!this.done) {
        this.done = true;
        this.array = undefined;
        this.ret.done = true;
        this.ret.value = undefined;
        arrayiters.cur--;
      }

      return this.ret;
    }

    return() {
      return this.finish();
    }
  }

  for (let i=0; i<arrayiters.length; i++) {
    arrayiters[i] = new ArrayIter();
  }

  Array.prototype[Symbol.iterator] = function() {
    return arrayiters[arrayiters.cur++].reset(this);
  }
}
