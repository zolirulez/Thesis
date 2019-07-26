function hBPf = hBPf_simulink(DQcor,DmQ,CRIT,hHR,dG)

DmIT = dG*1e-4*CRIT*48;
DmV = DmQ + DmIT;
hBPf = hHR-DQcor/DmV;

