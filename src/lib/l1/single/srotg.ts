/*

  jack dongarra, linpack, 3/11/78.
  jacob bogers, js port, 03/2018
*/

const { abs, sqrt } = Math;

export function srotg(p: { sa: number, sb: number, c: number, s: number }): void {

      let sa = p.sa;
      let sb = p.sb;
      const absA = abs(sa);
      const absB = abs(sb);

      let Z = 0;


      const roe = absA > absB ? sa : sb;
      let scale = absA + absB;

      if (scale === 0.0) {
            p.c = 1.0;
            p.s = 0.0;
            p.sa = 0; //R
            p.sb = 0; //Z
            return;
      }
      const v1 = sa / scale;
      const v2 = sb / scale;
      // scale > 0 so R > 0
      let R = scale * sqrt(v1 * v1 + v2 * v2);
      //R = sign(1.0, roe) * R;
      R = roe >= 0 ? R : -R;

      p.c = sa / R;
      p.s = sb / R;
      Z = 1.0;
      if (absA > absB) Z = p.s;
      if (absA <= absB && p.c !== 0.0) Z = 1.0 / p.c;
      // keep as is below
      // end original source code
      // repackage for exit
      p.sa = R;
      p.sb = Z;
}
