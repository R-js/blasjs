/* This is a conversion from BLAS to Typescript/Javascript
Copyright (C) 2018  Jacob K.F. Bogers  info@mail.jacob-bogers.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

import {
    Complex,
    errMissingIm,
    errWrongArg,
    FortranArr,
    isZero,
    lowerChar
} from '../../f_func';



export function chpr2(
    uplo: 'u' | 'l',
    n: number,
    alpha: Complex,
    x: FortranArr,
    incx: number,
    y: FortranArr,
    incy: number,
    ap: FortranArr
): void {

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }

    if (ap.i === undefined) {
        throw new Error(errMissingIm('ap.i'));
    }

    // faster then String.toLowerCase()
    const ul = lowerChar(uplo);

    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 5;
    }
    else if (incy === 0) {
        info = 7;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('chpr2', info));
    }
    const { re: AlphaRe, im: AlphaIm } = alpha;
    const alphaIsZero = isZero(alpha);

    if (n === 0 || alphaIsZero) return;

    let kx = incx < 0 ? 1 - (n - 1) * incx : 1;
    let ky = incy < 0 ? 1 - (n - 1) * incy : 1;

    let jx = kx;
    let jy = ky;

    let kk = 1;

    if (ul === 'u') {
        for (let j = 1; j <= n; j++) {
            //const apkk = kk - ap.base; //note kk is updated within this loop
            const xIsZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
            const yIsZero = y.r[jy - x.base] === 0 && y.i[jy - y.base] === 0;
            if (!xIsZero || !yIsZero) {
                //(a+ib)(c-id)=(ac+bd)+i(-ad+ bc)
                //TEMP1 = ALPHA*CONJG(Y(JY))

                let temp1Re = AlphaRe * y.r[jy - y.base] + AlphaIm * y.i[jy - y.base];
                let temp1Im = -AlphaRe * y.i[jy - y.base] + AlphaIm * y.r[jy - y.base];
                //(a-ib)(c+id)=(ac+bd)+i(ad-bc)
                //TEMP2 = CONJG(ALPHA*X(JX))

                let temp2Re = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
                let temp2Im = -(AlphaRe * x.i[jx - x.base] + AlphaIm * x.r[jx - x.base]);
                //console.log(`j${j}, t1=(${temp1Re},${temp1Im}), t2=(${temp2Re},${temp2Im})`);
                let ix = kx;
                let iy = ky;
                for (let k = kk; k <= kk + j - 2; k++) {
                    //AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                    //X(IX)*TEMP1 
                    // const xre = x.r[ix - x.base];
                    // const xim = x.i[ix - x.base];
                    // console.log(`j:${j},k:${k},ix:${ix},iy:${iy},x:(${xre},${xim}),`);

                    let re = (x.r[ix - x.base] * temp1Re - x.i[ix - x.base] * temp1Im);
                    let im = (x.r[ix - x.base] * temp1Im + x.i[ix - x.base] * temp1Re);
                    // Y(IY)*TEMP2
                    re += (y.r[iy - y.base] * temp2Re - y.i[iy - y.base] * temp2Im);
                    im += (y.r[iy - y.base] * temp2Im + y.i[iy - y.base] * temp2Re);

                    ap.r[k - ap.base] += re;
                    ap.i[k - ap.base] += im;

                    ix += incx;
                    iy += incy;
                }
                ap.i[kk + j - 1 - ap.base] = 0;
                ap.r[kk + j - 1 - ap.base] +=
                    //  REAL(X(JX)*TEMP1+Y(JY)*TEMP2)   
                    (x.r[jx - x.base] * temp1Re - x.i[jx - x.base] * temp1Im)
                    + //REAL(y[jy-y.base]*temp2)
                    (y.r[jy - y.base] * temp2Re - y.i[jy - y.base] * temp2Im);

            } //if  !xIsZero || !yIsZero
            else {
                // nuke imaginary part
                ap.i[kk + j - 1 - ap.base] = 0;
            }// if
            jx += incx;
            jy += incy;
            kk += j;
        }//for
    }
    else {
        // Form  A  when lower triangle is stored in AP.
        for (let j = 1; j <= n; j++) {
            //const apkk = kk - ap.base; //note kk is updated within this loop
            const xIsZero = x.r[jx - x.base] === 0 && x.i[jx - x.base] === 0;
            const yIsZero = y.r[jy - y.base] === 0 && y.i[jy - y.base] === 0;

            if (!(xIsZero && yIsZero)) {
                //TEMP1 = ALPHA*CONJG(Y(JY))
                //(a+ib)(c-id)=(ac+bd)+i(-ad+ bc)
                let temp1Re = AlphaRe * y.r[jy - y.base] + AlphaIm * y.i[jy - y.base];
                let temp1Im = -AlphaRe * y.i[jy - y.base] + AlphaIm * y.r[jy - y.base];

                //(a-ib)(c+id)=(ac+bd)+i(ad-bc)
                //TEMP2 = CONJG(ALPHA*X(JX))
                let temp2Re = AlphaRe * x.r[jx - x.base] - AlphaIm * x.i[jx - x.base];
                let temp2Im = -(AlphaRe * x.i[jx - x.base] + AlphaIm * x.r[jx - x.base]);

                ap.r[kk - ap.base] += (x.r[jx - x.base] * temp1Re - x.i[jx - x.base] * temp1Im)
                    + (y.r[jy - y.base] * temp2Re - y.i[jy - y.base] * temp2Im);


                ap.i[kk - ap.base] = 0;

                let ix = jx;
                let iy = jy;

                for (let k = kk + 1; k <= kk + n - j; k++) {

                    ix += incx;
                    iy += incy;
                    // AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2

                    ap.r[k - ap.base] +=
                        (x.r[ix - x.base] * temp1Re - x.i[ix - x.base] * temp1Im)
                        +
                        (y.r[iy - y.base] * temp2Re - y.i[iy - y.base] * temp2Im);
                    ap.i[k - ap.base] +=
                        (x.r[ix - x.base] * temp1Im + x.i[ix - x.base] * temp1Re)
                        +
                        (y.r[iy - y.base] * temp2Im + y.i[iy - y.base] * temp2Re);
                }
            } else {
                ap.i[kk - ap.base] = 0;
            }
            jx += incx;
            jy += incy;
            kk += n - j + 1;
        }
    }
}

