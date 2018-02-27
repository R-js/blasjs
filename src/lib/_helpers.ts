/*
Intrinsic Fortran functions
Sign(A, B)
Sign: INTEGER or REAL function, the exact type being the result of cross - promoting the types of all the arguments.

    A: INTEGER or REAL; scalar; INTENT(IN).

        B: INTEGER or REAL; scalar; INTENT(IN).

Intrinsic groups: (standard FORTRAN 77).

Description:

Returns`ABS(A)*s', where s is +1 if `B.GE.0', -1 otherwise.

See Abs Intrinsic, for the function that returns the magnitude of a value.
*/

export function fsign(a: number, b?: number): number {
    if (b === undefined) {
        return a;
    }
    const rc = Math.abs(a);
    return b >= 0 ? rc : -rc;
}
