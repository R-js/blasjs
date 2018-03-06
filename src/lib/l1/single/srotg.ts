/*

  jack dongarra, linpack, 3/11/78.
  jacob bogers, js port, 03/2018
*/

const { abs, sqrt, pow } = Math;
import { sign } from '../../f_func';

export function srotg(p: { sa: number, sb: number, c: number, s: number }): void {

      let SA = p.sa;
      let SB = p.sb;
      let C = p.c;
      let S = p.s;

      let roe = SB;
      let Z = 0;
      let R: number;

      if (abs(SA) > abs(SB)) roe = SA;
      let scale = abs(SA) + abs(SB);
      if (scale === 0.0) {
            C = 1.0;
            S = 0.0;
            R = 0.0;
            Z = 0.0;
      }
      else {
            R = scale * sqrt(pow(SA / scale, 2) + pow(SB / scale, 2));
            R = sign(1.0, roe) * R;
            C = SA / R;
            S = SB / R;
            Z = 1.0;
            if (abs(SA) > abs(SB)) Z = S;
            if (abs(SB) >= abs(SA) && C !== 0.0) Z = 1.0 / C
      }
      // keep as is below
      SA = R;
      SB = Z;
      // end original source code
      // repackage for exit
      p.sa = SA;
      p.sb = SB;
      p.c = C;
      p.s = S;
}
