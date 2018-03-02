//intrinsic routines of fortran

export function sign(a: number, b?: number): number {
    if (b === undefined) {
        return a;
    }
    const rc = Math.abs(a);
    return b >= 0 ? rc : -rc;
}
