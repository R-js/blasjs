   log('raw gpu version:');
    const canvas = document.createElement('canvas');
    document.body.appendChild(canvas);
    const webgl = canvas.getContext('webgl2');
    const gpuKernel = gpuFunctionRaw(webgl, canvas);
    log(gpuKernel(a, b, c, lda, ldb, ldc, m, n, k, alpha, beta));


const gpuFunctionRaw = (gl, canvas, debug) => {
  const paramNames = ["a","b","c","lda","ldb","ldc","m","n","k","alpha","beta"];
  const paramTypes = ["Array","Array","Array","Integer","Integer","Integer","Integer","Integer","Integer","Float","Float"];
  const texSize = [4,6];
  const output = [4,6];
  const outputDim = [4,6,1];
  const compiledFragShaderString = `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;
const float LOOP_MAX = 1000.0;
uniform ivec3 uOutputDim;
uniform ivec2 uTexSize;
in vec2 vTexCoord;
vec2 integerMod(vec2 x, float y) {
  vec2 res = floor(mod(x, y));
  return res * step(1.0 - floor(y), -res);
}
vec3 integerMod(vec3 x, float y) {
  vec3 res = floor(mod(x, y));
  return res * step(1.0 - floor(y), -res);
}
vec4 integerMod(vec4 x, vec4 y) {
  vec4 res = floor(mod(x, y));
  return res * step(1.0 - floor(y), -res);
}
float integerMod(float x, float y) {
  float res = floor(mod(x, y));
  return res * (res > floor(y) - 1.0 ? 0.0 : 1.0);
}
int integerMod(int x, int y) {
  return x - (y * int(x/y));
}
			  float div_with_int_check(float x, float y) {
			  if (floor(x) == x && floor(y) == y && integerMod(x, y) == 0.0) {
			    return float(int(x)/int(y));
			  }
			  return x / y;
			}
// Here be dragons!
// DO NOT OPTIMIZE THIS CODE
// YOU WILL BREAK SOMETHING ON SOMEBODY'S MACHINE
// LEAVE IT AS IT IS, LEST YOU WASTE YOUR OWN TIME
const vec2 MAGIC_VEC = vec2(1.0, -256.0);
const vec4 SCALE_FACTOR = vec4(1.0, 256.0, 65536.0, 0.0);
const vec4 SCALE_FACTOR_INV = vec4(1.0, 0.00390625, 0.0000152587890625, 0.0); // 1, 1/256, 1/65536
float decode32(vec4 rgba) {
  rgba *= 255.0;
  vec2 gte128;
  gte128.x = rgba.b >= 128.0 ? 1.0 : 0.0;
  gte128.y = rgba.a >= 128.0 ? 1.0 : 0.0;
  float exponent = 2.0 * rgba.a - 127.0 + dot(gte128, MAGIC_VEC);
  float res = exp2(round(exponent));
  rgba.b = rgba.b - 128.0 * gte128.x;
  res = dot(rgba, SCALE_FACTOR) * exp2(round(exponent-23.0)) + res;
  res *= gte128.y * -2.0 + 1.0;
  return res;
}
vec4 encode32(float f) {
  float F = abs(f);
  float sign = f < 0.0 ? 1.0 : 0.0;
  float exponent = floor(log2(F));
  float mantissa = (exp2(-exponent) * F);
  // exponent += floor(log2(mantissa));
  vec4 rgba = vec4(F * exp2(23.0-exponent)) * SCALE_FACTOR_INV;
  rgba.rg = integerMod(rgba.rg, 256.0);
  rgba.b = integerMod(rgba.b, 128.0);
  rgba.a = exponent*0.5 + 63.5;
  rgba.ba += vec2(integerMod(exponent+127.0, 2.0), sign) * 128.0;
  rgba = floor(rgba);
  rgba *= 0.003921569; // 1/255
  return rgba;
}
// Dragons end here
float decode(vec4 rgba, int x, int bitRatio) {
  if (bitRatio == 1) {
    return decode32(rgba);
  }
  int channel = integerMod(x, bitRatio);
  if (bitRatio == 4) {
    return rgba[channel] * 255.0;
  }
  else {
    return rgba[channel*2] * 255.0 + rgba[channel*2 + 1] * 65280.0;
  }
}
int index;
ivec3 threadId;
ivec3 indexTo3D(int idx, ivec3 texDim) {
  int z = int(idx / (texDim.x * texDim.y));
  idx -= z * int(texDim.x * texDim.y);
  int y = int(idx / texDim.x);
  int x = int(integerMod(idx, texDim.x));
  return ivec3(x, y, z);
}
float get(sampler2D tex, ivec2 texSize, ivec3 texDim, int bitRatio,  int z, int y, int x) {
  ivec3 xyz = ivec3(x, y, z);
  int index = xyz.x + texDim.x * (xyz.y + texDim.y * xyz.z);
  int w = texSize.x;
  vec2 st = vec2(float(integerMod(index, w)), float(index / w)) + 0.5;
  vec4 texel = texture(tex, st / vec2(texSize));
  return decode(texel, x, bitRatio);
}
vec4 getImage2D(sampler2D tex, ivec2 texSize, ivec3 texDim, int z, int y, int x) {
  ivec3 xyz = ivec3(x, y, z);
  int index = xyz.x + texDim.x * (xyz.y + texDim.y * xyz.z);
  int w = texSize.x;
  vec2 st = vec2(float(integerMod(index, w)), float(index / w)) + 0.5;
  return texture(tex, st / vec2(texSize));
}
vec4 getImage3D(sampler2DArray tex, ivec2 texSize, ivec3 texDim, int z, int y, int x) {
  ivec3 xyz = ivec3(x, y, z);
  int index = xyz.x + texDim.x * (xyz.y + texDim.y * xyz.z);
  int w = texSize.x;
  vec2 st = vec2(float(integerMod(index, w)), float(index / w)) + 0.5;
  return texture(tex, vec3(st / vec2(texSize), z));
}
float get(sampler2D tex, ivec2 texSize, ivec3 texDim, int bitRatio, int y, int x) {
  return get(tex, texSize, texDim, bitRatio, 0, y, x);
}
float get(sampler2D tex, ivec2 texSize, ivec3 texDim, int bitRatio, int x) {
  return get(tex, texSize, texDim, bitRatio, 0, 0, x);
}
vec4 getImage2D(sampler2D tex, ivec2 texSize, ivec3 texDim, int y, int x) {
  return getImage2D(tex, texSize, texDim, 0, y, x);
}
vec4 getImage2D(sampler2D tex, ivec2 texSize, ivec3 texDim, int x) {
  return getImage2D(tex, texSize, texDim, 0, 0, x);
}
vec4 actualColor;
void color(float r, float g, float b, float a) {
  actualColor = vec4(r,g,b,a);
}
void color(float r, float g, float b) {
  color(r,g,b,1.0);
}
uniform highp sampler2D user_a;
uniform highp ivec2 user_aSize;
uniform highp ivec3 user_aDim;
uniform highp int user_aBitRatio;
uniform highp sampler2D user_b;
uniform highp ivec2 user_bSize;
uniform highp ivec3 user_bDim;
uniform highp int user_bBitRatio;
uniform highp sampler2D user_c;
uniform highp ivec2 user_cSize;
uniform highp ivec3 user_cDim;
uniform highp int user_cBitRatio;
uniform float user_lda;
uniform float user_ldb;
uniform float user_ldc;
uniform float user_m;
uniform float user_n;
uniform float user_k;
uniform float user_alpha;
uniform float user_beta;
out vec4 data0;
float kernelResult = 0.0;
void kernel() {
float user_cIndex=((user_ldc*float(threadId.y))+float(threadId.x));
float user_sum=(get(user_c, ivec2(user_cSize[0],user_cSize[1]), ivec3(user_cDim[0],user_cDim[1],user_cDim[2]), user_cBitRatio, int(user_cIndex))*user_beta);
for (float user_i=0.0;user_i<LOOP_MAX;user_i++){
if (user_i<user_k) {
float user_aIndex=((user_lda*user_i)+float(threadId.x));float user_bIndex=((user_ldb*float(threadId.y))+user_i);user_sum+=((get(user_a, ivec2(user_aSize[0],user_aSize[1]), ivec3(user_aDim[0],user_aDim[1],user_aDim[2]), user_aBitRatio, int(user_aIndex))*user_alpha)*get(user_b, ivec2(user_bSize[0],user_bSize[1]), ivec3(user_bDim[0],user_bDim[1],user_bDim[2]), user_bBitRatio, int(user_bIndex)));
} else {
break;
}
}
kernelResult = user_sum;return;
}
void main(void) {
  index = 
     int(vTexCoord.s * float(uTexSize.x)) + 
     int(vTexCoord.t * float(uTexSize.y)) * uTexSize.x;
  threadId = indexTo3D(index, uOutputDim);
  kernel();
  data0 = encode32(kernelResult);
}`;
  const compiledVertShaderString = `#version 300 es
precision highp float;
precision highp int;
precision highp sampler2D;
in vec2 aPos;
in vec2 aTexCoord;
out vec2 vTexCoord;
uniform vec2 ratio;
void main(void) {
  gl_Position = vec4((aPos + vec2(1)) * ratio + vec2(-1), 0, 1);
  vTexCoord = aTexCoord;
}`;
  const programUniformLocationCache = {};
  const textureCache = {};
  const uniform1fCache = {};
  const uniform1iCache = {};
  const uniform2fCache = {};
  const uniform2ivCache = {};
  const uniform3ivCache = {};

  // setup functions
  function clone(obj) {
    if (obj === null || (typeof obj === 'undefined' ? 'undefined' : typeof(obj)) !== 'object' || obj.hasOwnProperty('isActiveClone')) return obj;

    var temp = obj.constructor();

    for (var key in obj) {
      if (Object.prototype.hasOwnProperty.call(obj, key)) {
        obj.isActiveClone = null;
        temp[key] = clone(obj[key]);
        delete obj.isActiveClone;
      }
    }
    return temp;
  }
  function isArray(array) {
    if (isNaN(array.length)) {
      return false;
    }
    return true;
  }
  function flattenTo(array, target) {
    if (isArray(array[0])) {
      if (isArray(array[0][0])) {
        flatten3dArrayTo(array, target);
      } else {
        flatten2dArrayTo(array, target);
      }
    } else {
      target.set(array);
    }
  }
  function flatten2dArrayTo(array, target) {
    var offset = 0;
    for (var y = 0; y < array.length; y++) {
      target.set(array[y], offset);
      offset += array[y].length;
    }
  }
  function flatten3dArrayTo(array, target) {
    var offset = 0;
    for (var z = 0; z < array.length; z++) {
      for (var y = 0; y < array[z].length; y++) {
        target.set(array[z][y], offset);
        offset += array[z][y].length;
      }
    }
  }
  function formatArrayTransfer(value, length) {
    var bitRatio = 1;
    var valuesFlat = new Float32Array(length);
    flattenTo(value, valuesFlat);
    return {
      bitRatio: bitRatio,
      valuesFlat: valuesFlat
    };
  }
  function getDimensions(x, pad) {
    var ret = void 0;
    if (isArray(x)) {
      var dim = [];
      var temp = x;
      while (isArray(temp)) {
        dim.push(temp.length);
        temp = temp[0];
      }
      ret = dim.reverse();
    } else {
      throw 'Unknown dimensions of ' + x;
    }

    if (pad) {
      ret = clone(ret);
      while (ret.length < 3) {
        ret.push(1);
      }
    }
    return new Int32Array(ret);
  }
  function dimToTexSize(opt, dimensions, output) {
    var numTexels = dimensions[0];
    var w = dimensions[0];
    var h = dimensions[1];
    for (var i = 1; i < dimensions.length; i++) {
      numTexels *= dimensions[i];
    }

    if (opt.floatTextures && (!output || opt.floatOutput)) {
      w = numTexels = Math.ceil(numTexels / 4);
    }
    if (h > 1 && w * h === numTexels) {
      return [w, h];
    }
    var sqrt = Math.sqrt(numTexels);
    var high = Math.ceil(sqrt);
    var low = Math.floor(sqrt);
    while (high * low > numTexels) {
      high--;
      low = Math.ceil(numTexels / high);
    }
    w = low;
    h = Math.ceil(numTexels / w);
    return [w, h];
  }
  function setUniform2f(name, value1, value2) {
    if (uniform2fCache.hasOwnProperty(name)) {
      var cache = uniform2fCache[name];
      if (value1 === cache[0] && value2 === cache[1]) {
        return;
      }
    }
    uniform2fCache[name] = [value1, value2];
    var loc = getUniformLocation(name);
    gl.uniform2f(loc, value1, value2);
  }
  function setUniform1i(name, value) {
    if (uniform1iCache.hasOwnProperty(name)) {
      var cache = uniform1iCache[name];
      if (value === cache) {
        return;
      }
    }
    uniform1iCache[name] = value;
    var loc = getUniformLocation(name);
    gl.uniform1i(loc, value);
  }
  function setUniform1f(name, value) {
    if (uniform1fCache.hasOwnProperty(name)) {
      var cache = uniform1fCache[name];
      if (value === cache) {
        return;
      }
    }
    uniform1fCache[name] = value;
    var loc = getUniformLocation(name);
    gl.uniform1f(loc, value);
  }
  function setUniform2iv(name, value) {
    if (uniform2ivCache.hasOwnProperty(name)) {
      var cache = uniform2ivCache[name];
      if (value[0] === cache[0] && value[1] === cache[1]) {
        return;
      }
    }
    uniform2ivCache[name] = value;
    var loc = getUniformLocation(name);
    gl.uniform2iv(loc, value);
  }
  function setUniform3iv(name, value) {
    if (uniform3ivCache.hasOwnProperty(name)) {
      var cache = uniform3ivCache[name];
      if (value[0] === cache[0] && value[1] === cache[1] && value[2] === cache[2]) {
        return;
      }
    }
    uniform3ivCache[name] = value;
    var loc = getUniformLocation(name);
    gl.uniform3iv(loc, value);
  }
  function getUniformLocation(name) {
    if (programUniformLocationCache.hasOwnProperty(name)) {
      return programUniformLocationCache[name];
    }
    return programUniformLocationCache[name] = gl.getUniformLocation(program, name);
  }
  function getTexture(name) {
    if (textureCache.hasOwnProperty(name)) {
      return textureCache[name];
    }
    return textureCache[name] = gl.createTexture();
  }
  function splitArray(array, part) {
    var result = [];
    for (var i = 0; i < array.length; i += part) {
      result.push(new array.constructor(array.buffer, i * 4 + array.byteOffset, part));
    }
    return result;
  }
  gl.enable(gl.SCISSOR_TEST);
  gl.viewport(0, 0, texSize[0], texSize[1]);
  canvas.width = texSize[0];
  canvas.height = texSize[1];

  var vertShader = gl.createShader(gl.VERTEX_SHADER);
  gl.shaderSource(vertShader, compiledVertShaderString);
  gl.compileShader(vertShader);

  var fragShader = gl.createShader(gl.FRAGMENT_SHADER);
  gl.shaderSource(fragShader, compiledFragShaderString);
  gl.compileShader(fragShader);

  if (!gl.getShaderParameter(vertShader, gl.COMPILE_STATUS)) {
    console.log(compiledVertShaderString);
    console.error('An error occurred compiling the shaders: ' + gl.getShaderInfoLog(vertShader));
    throw new Error('Error compiling vertex shader');
  }
  if (!gl.getShaderParameter(fragShader, gl.COMPILE_STATUS)) {
    console.log(compiledFragShaderString);
    console.error('An error occurred compiling the shaders: ' + gl.getShaderInfoLog(fragShader));
    throw new Error('Error compiling fragment shader');
  }

  if (debug) {
    console.log('Options:');
    console.log('GLSL Shader Output:');
    console.log(compiledFragShaderString);
  }

  var program = gl.createProgram();
  gl.attachShader(program, vertShader);
  gl.attachShader(program, fragShader);
  gl.linkProgram(program);
  var framebuffer = gl.createFramebuffer();
  framebuffer.width = texSize[0];
  framebuffer.height = texSize[1];

  var vertices = new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]);
  var texCoords = new Float32Array([0, 0, 1, 0, 0, 1, 1, 1]);
  var texCoordOffset = vertices.byteLength;
  var buffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
  gl.bufferData(gl.ARRAY_BUFFER, vertices.byteLength + texCoords.byteLength, gl.STATIC_DRAW);
  gl.bufferSubData(gl.ARRAY_BUFFER, 0, vertices);
  gl.bufferSubData(gl.ARRAY_BUFFER, texCoordOffset, texCoords);

  var aPosLoc = gl.getAttribLocation(program, 'aPos');
  gl.enableVertexAttribArray(aPosLoc);
  gl.vertexAttribPointer(aPosLoc, 2, gl.FLOAT, gl.FALSE, 0, 0);

  var aTexCoordLoc = gl.getAttribLocation(program, 'aTexCoord');
  gl.enableVertexAttribArray(aTexCoordLoc);
  gl.vertexAttribPointer(aTexCoordLoc, 2, gl.FLOAT, gl.FALSE, 0, texCoordOffset);
  gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffer);

  //setup output
  var outputTexture = gl.createTexture();
  gl.activeTexture(gl.TEXTURE0 + paramNames.length);
  gl.bindTexture(gl.TEXTURE_2D, outputTexture);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, texSize[0], texSize[1], 0, gl.RGBA, gl.UNSIGNED_BYTE, null);
  gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, outputTexture, 0);

  return function() {
    gl.useProgram(program);
    gl.scissor(0, 0, texSize[0], texSize[1]);

    setUniform3iv('uOutputDim', new Int32Array(outputDim));
    setUniform2iv('uTexSize', texSize);
    setUniform2f('ratio', 1, 1);

    // setup arguments
    for (var texIndex = 0; texIndex < paramNames.length; texIndex++) {
      var value = arguments[texIndex];
      var type = paramTypes[texIndex];
      var name = paramNames[texIndex];
      switch (type) {
        case 'Array':
          var dim = getDimensions(value, true);
          var size = dimToTexSize({}, dim);
          gl.activeTexture(gl.TEXTURE0 + texIndex);
          gl.bindTexture(gl.TEXTURE_2D, getTexture('ARGUMENT_' + name));
          gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
          gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
          gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
          gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);

          var length = size[0] * size[1];
          var _formatArrayTransfer = formatArrayTransfer(value, length);
          var valuesFlat = _formatArrayTransfer.valuesFlat;
          var bitRatio = _formatArrayTransfer.bitRatio;
          var buffer = new Uint8Array(valuesFlat.buffer);
          gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, size[0] / bitRatio, size[1], 0, gl.RGBA, gl.UNSIGNED_BYTE, buffer);

          setUniform3iv('user_' + name + 'Dim', dim);
          setUniform2iv('user_' + name + 'Size', size);
          setUniform1i('user_' + name + 'BitRatio', bitRatio);
          setUniform1i('user_' + name, texIndex);
          break;
        case 'Integer':
        case 'Float':
          setUniform1f('user_' + name, value);
          break;
        default:
          throw new Error('unsupported type argument type ' + type);
      }
    }

    gl.bindFramebuffer(gl.FRAMEBUFFER, framebuffer);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);

    var bytes = new Uint8Array(texSize[0] * texSize[1] * 4);
    gl.readPixels(0, 0, texSize[0], texSize[1], gl.RGBA, gl.UNSIGNED_BYTE, bytes);
    var result = new Float32Array(bytes.buffer);
    result = result.subarray(0, outputDim[0] * outputDim[1] * outputDim[2]);

    if (output.length === 1) {
      return result;
    } else if (output.length === 2) {
      return splitArray(result, output[0]);
    } else if (output.length === 3) {
      var cube = splitArray(result, output[0] * output[1]);
      return cube.map(function (x) {
        return splitArray(x, output[0]);
      });
    }
    return result;
  };
};
