#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _Cad_reg(void);
extern void _CaH_reg(void);
extern void _CaT_reg(void);
extern void _charge_reg(void);
extern void _h_reg(void);
extern void _Kca_reg(void);
extern void _Kv1_axonal_reg(void);
extern void _Kv7_reg(void);
extern void _Kv_reg(void);
extern void _na8st_reg(void);
extern void _nax8st_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," HallermannEtAl2012/Cad.mod");
    fprintf(stderr," HallermannEtAl2012/CaH.mod");
    fprintf(stderr," HallermannEtAl2012/CaT.mod");
    fprintf(stderr," HallermannEtAl2012/charge.mod");
    fprintf(stderr," HallermannEtAl2012/h.mod");
    fprintf(stderr," HallermannEtAl2012/Kca.mod");
    fprintf(stderr," HallermannEtAl2012/Kv1_axonal.mod");
    fprintf(stderr," HallermannEtAl2012/Kv7.mod");
    fprintf(stderr," HallermannEtAl2012/Kv.mod");
    fprintf(stderr," HallermannEtAl2012/na8st.mod");
    fprintf(stderr," HallermannEtAl2012/nax8st.mod");
    fprintf(stderr, "\n");
  }
  _Cad_reg();
  _CaH_reg();
  _CaT_reg();
  _charge_reg();
  _h_reg();
  _Kca_reg();
  _Kv1_axonal_reg();
  _Kv7_reg();
  _Kv_reg();
  _na8st_reg();
  _nax8st_reg();
}
