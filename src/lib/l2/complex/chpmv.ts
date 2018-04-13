import {
    Complex,
    errMissingIm,
    errWrongArg,
    FortranArr,
    isOne,
    isZero,
    lowerChar,
    mul_cxc,
    mul_cxr,
    mul_rxc,
    mul_rxr
} from '../../f_func';

/**
*>  -- Jacob Bogers, 2018/03, jkfbogers@gmail.com
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs. 
*/
/* 
    y := alpha*A*x + beta*y,
    A is an nxn matrix (upper or lower) in packed form
*/
export function chpmv(
    uplo: 'u' | 'l',
    n: number,
    alpha: Complex,
    ap: FortranArr,
    x: FortranArr,
    incx: number,
    beta: Complex,
    y: FortranArr,
    incy: number
): void {

    if (x.i === undefined) {
        throw new Error(errMissingIm('x.i'));
    }

    if (y.i === undefined) {
        throw new Error(errMissingIm('y.i'));
    }

    if (ap.i === undefined) {
        throw new Error(errMissingIm('a.i'));
    }

    const ul = lowerChar(uplo);


    let info = 0;
    if (!'ul'.includes(ul)) {
        info = 1;
    }
    else if (n < 0) {
        info = 2;
    }
    else if (incx === 0) {
        info = 3;
    }
    else if (incy === 0) {
        info = 9;
    }
    if (info !== 0) {
        throw new Error(errWrongArg('chpmv', info));
    }


    const betaIsOne = isOne(beta);
    const alphaIsZero = isZero(alpha);
    const betaIsZero = isZero(beta);



    if (n === 0 || (alphaIsZero && betaIsOne)) return;

    let kx = incx > 0 ? 1 : 1 - (n - 1) * incx;
    let ky = incy > 0 ? 1 : 1 - (n - 1) * incy;

    // First form  y := beta*y.
    let iy = ky;
    //console.log({ betaIsOne });
    if (!betaIsOne) {
        for (let i = 1; i <= n; i++) {
            const re = betaIsZero ? 0 : beta.re * y.r[iy - y.base] - beta.im * y.i[iy - y.base];
            const im = betaIsZero ? 0 : beta.re * y.i[iy - y.base] + beta.im * y.r[iy - y.base];
            y.r[iy - y.base] = re;
            y.i[iy - y.base] = im;
            iy += incy;
        }
    }
    // }
    //y is ok here
    //(y.toArr() as Complex[]).forEach(c => console.log(`(${c.re},${c.im})`));
    if (alphaIsZero) return;
    let kk = 1;
    if (ul === 'u') {
        //Form  y  when AP contains the upper triangle.
        let jx = kx;
        let jy = ky;

        for (let j = 1; j <= n; j++) {
            let temp1Re = alpha.re * x.r[jx - x.base] - alpha.im * x.i[jx - x.base];
            let temp1Im = alpha.re * x.i[jx - x.base] + alpha.im * x.r[jx - x.base];
            //console.log(`${j}, (${temp1Re},${temp1Im})`);
            let temp2Re = 0;
            let temp2Im = 0;

            let ix = kx;
            let iy = ky;

            for (let k = kk; k <= kk + j - 2; k++) {
                const apk = k - ap.base;
                // Y(IY) = Y(IY) + TEMP1*AP(K)
                y.r[iy - y.base] += temp1Re * ap.r[apk] - temp1Im * ap.i[apk];
                y.i[iy - y.base] += temp1Re * ap.i[apk] + temp1Im * ap.r[apk];
                // TEMP2 = TEMP2 + DCONJG(AP(K))*X(IX)
                // (a-ib)*(c+id) = (ac+bd)+i(ad-bc)
                temp2Re += ap.r[apk] * x.r[ix - x.base] + ap.i[apk] * x.i[ix - x.base];
                temp2Im += ap.r[apk] * x.i[ix - x.base] - ap.i[apk] * x.r[ix - x.base];
                //  console.log(`${ix},${iy}, (${y.r[iy - y.base]},${y.i[iy - y.base]})`);
                ix += incx;
                iy += incy;
            }
            // Y(JY) + TEMP1*DBLE(AP(KK+J-1)) + ALPHA*TEMP2
            const dbleapkj = ap.r[kk + j - 1 - ap.base];
            const re1 = temp1Re * dbleapkj;
            const im1 = temp1Im * dbleapkj;

            const re2 = (alpha.re * temp2Re - alpha.im * temp2Im);
            const im2 = (alpha.re * temp2Im + alpha.im * temp2Re);


            y.r[jy - y.base] += re1 + re2;
            y.i[jy - y.base] += im1 + im2;
            //console.log(`${jy},(${y.r[jy - y.base]},${y.i[jy - y.base]}) `);
            jx += incx;
            jy += incy;
            kk += j;
        }
    }
    else {
        //*        Form  y  when AP contains the lower triangle.
        let jx = kx;
        let jy = ky;
        for (let j = 1; j <= n; j++) {
            // TEMP1 = ALPHA*X(JX)
            const c0 = mul_cxr(
                alpha,
                x.r[jx - x.base],
                x.i[jx - x.base]);
            let temp1Re = c0.re;
            let temp1Im = c0.im;

            let temp2Re = 0;
            let temp2Im = 0;
            const dbleapkk = ap.r[kk - ap.base]
            //   Y(JY) = Y(JY) + TEMP1*DBLE(AP(KK))
            y.r[jy - y.base] += temp1Re * dbleapkk;
            y.i[jy - y.base] += temp1Im * dbleapkk;
            //            console.log(`${jy} ( ${y.r[jy - y.base]},${y.i[jy - y.base]}`)
            let ix = jx;
            let iy = jy;
            for (let k = kk + 1; k <= kk + n - j; k++) {
                ix += incx;
                iy += incy;
                const apk = k - ap.base;
                //    Y(IY) = Y(IY) + TEMP1*AP(K)
                //console.log(`${k}, ${ix}, ${iy}, (${y.r[iy - y.base]},${y.r[iy - y.base]}) `);
                const c1 = mul_rxr(temp1Re, temp1Im, ap.r[apk], ap.i[apk]);
                y.r[iy - y.base] += c1.re;
                y.i[iy - y.base] += c1.im;
                //console.log(`${k}, ${ix}, ${iy}, (${y.r[iy - y.base]},${y.r[iy - y.base]}) `);

                //  TEMP2 = TEMP2 + DCONJG(AP(K))*X(IX)
                // (a-ib)(c+id) = (ac+bd)+i(ad-bc)
                const c2 = mul_rxr(
                    ap.r[apk], -ap.i[apk],
                    x.r[ix - x.base], x.i[ix - x.base]);
                temp2Re += c2.re;
                temp2Im += c2.im;
            }
            // Y(JY) = Y(JY) + ALPHA*TEMP2
            const c3 = mul_cxr(
                alpha,
                temp2Re, temp2Im
            );
            y.r[jy - y.base] += c3.re;
            y.i[jy - y.base] += c3.im;
            jx += incx;
            jy += incy;
            kk += (n - j + 1);
        }
    }
}

