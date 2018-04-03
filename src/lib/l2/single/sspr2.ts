/** 
 -- Jacob Bogers, 03/2018, jkfbogers@gmail.com
 
 -- Written on 22-October-1986.
     Jack Dongarra, Argonne National Lab.
     Jeremy Du Croz, Nag Central Office.
     Sven Hammarling, Nag Central Office.
     Richard Hanson, Sandia National Labs.
*/

import { errWrongArg, FortranArr, lowerChar } from '../../f_func';

/*
 SSYMV  performs the matrix-vector  operation

 
    
    A := alpha*x*y**T + alpha*y*x**T + A,
    
*/

export function sspr2(
    uplo: 'u' | 'l',
    n: number,
    alpha: number,
    x: FortranArr,
    incx: number,
    y: FortranArr,
    incy: number,
    ap: FortranArr): void {

    // validate input parameters    

    let info = 0;
    const ul = lowerChar(uplo);
    //console.log('start', ul);

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
        throw new Error(errWrongArg('sspr2', info));
    }
    //console.log('startw2');
    //     Quick return if possible.

    if (n === 0 || alpha === 0) return;

    let kx = (incx > 0) ? 1 : 1 - (n - 1) * incx;
    let ky = (incy > 0) ? 1 : 1 - (n - 1) * incy;

    let jx = kx;
    let jy = ky;
    //console.log('startw3');
    //    Start the operations. In this version the elements of the array AP
    //    are accessed sequentially with one pass through AP.

    let kk = 1;
    if (ul === 'u') {
        //Form  A  when upper triangle is stored in AP.
        for (let j = 1; j <= n; j++) {
            if (x.r[jx - x.base] !== 0 || y.r[jy - y.base] !== 0) {
                let temp1 = alpha * y.r[jy - y.base];
                //console.log({ temp1 });
                let temp2 = alpha * x.r[jx - x.base];
                //console.log({ temp2 });
                let ix = kx;
                let iy = ky;
                for (let k = kk; k <= kk + j - 1; k++) {
                    ap.r[k - ap.base] += x.r[ix - x.base] * temp1 + y.r[iy - y.base] * temp2;
                    ix += incx;
                    iy += incy;
                }
            }
            jx += incx;
            jy += incy;
            kk += j;
        }
    }
    else {
        // Form  A  when lower triangle is stored in AP.
        console.log('startw4');
        for (let j = 1; j <= n; j++) {
            //console.log('startw5');
            const xAndyIsZero = x.r[jx - x.base] === 0 && y.r[jy - y.base] === 0;
            //console.log(`xyIsZero=${xAndyIsZero}`);
            if (!xAndyIsZero) {
                //console.log('startw6', alpha, y.r);
                let temp1 = alpha * y.r[jy - y.base];
                let temp2 = alpha * x.r[jx - x.base];
                //console.log({ temp1 });
                let ix = jx;
                let iy = jy;
                for (let k = kk; k <= kk + n - j; k++) {
                    ap.r[k - ap.base] += x.r[ix - x.base] * temp1 + y.r[iy - y.base] * temp2;
                    //console.log({ temp1 });
                    ix += incx;
                    iy += incy;
                }
            }
            //console.log('startw7');
            jx += incx;
            jy += incy;
            kk += n - j + 1;
        }
    }
}
