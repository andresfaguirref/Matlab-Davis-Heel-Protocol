function punto = L_RC_LocalValue(dis_asis, ldepth, leg_length)

punto1 = [ldepth, dis_asis, ldepth];

punto2 = punto1.*[-0.951, 0.5, -0.309];

%aparte

punto3 = leg_length*0.115;
punto3 = punto3-0.015;

punto4  =punto3.*[0.272, -0.476, -0.837];

punto = punto4+punto2;

end