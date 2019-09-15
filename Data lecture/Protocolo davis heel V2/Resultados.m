
ni=316;
nf=709;

Xx = Find_name(ProcesadosBTS_names, 'LA.X');
Yy = Find_name(ProcesadosBTS_names, 'LA.Y');
Zz = Find_name(ProcesadosBTS_names, 'LA.Z');


Gold_LA_RC = [ProcesadosBTS_3Dpoints(ni:nf,Xx), ProcesadosBTS_3Dpoints(ni:nf,Yy), ProcesadosBTS_3Dpoints(ni:nf,Zz)];

resta_LA = Gold_LA_RC - RC_left_ankle;