
tab = readtable('Gasloop_20190606.Csv');
addpath('C:\Users\u375749\Documents\Thesis\Codes\ModelicaPlant')
tab.Properties.VariableNames'
figure(1)
start = 1.5e4;
plot(tab{start:end,2:end})
fielddata = tab{start:end,:};
%     {'t'        }
%     {'pGC'      }
%     {'TA0'      }
%     {'pGCr'     }
%     {'TGC'      }
%     {'TGCr'     }
%     {'CRA'      }
%     {'CRV'      }
%     {'THR_wrong'}
%     {'TMT'      }
%     {'TMTi'     }
%     {'CRMT'     }
%     {'ToMT'     }
%     {'ToMTr'    }
%     {'pR'       }
%     {'TR'       }
%     {'pRr'      }
%     {'TRr'      }
%     {'CRG'      }
%     {'ToITr'    }
%     {'CRIT'     }
%     {'TIT'      }
%     {'TITi'     }
%     {'xiR'      }
%     {'ToLT'     }
%     {'ToLTr'    }
%     {'CRLT'     }
%     {'TLT'      }
%     {'TLTi'     }
% Parameters
VR = 0.115;
VMT = (6.5 + 9.2)/3600/48*2; % Assuming that maximal f is double of nominal 24 Hz.
VIT = 3.3/3600/48;
Mxf = 48;
delay = 75;
TA0offsetindex = 12365; % When the TA0 is correct again (else 5 minus is considered to find gasloop potential references.
% Data
CRA = fielddata(:,7)/100;
TBP = fielddata(:,5)+273.15;
TMT = fielddata(:,10)+273.15;
TMTd = tab{start-delay:end-delay,10}+273.15;
BP = zeros(length(fielddata),1);
CRV = fielddata(:,8)/100;
CRIT = fielddata(:,21)/100;
CRG = fielddata(:,19)/100;
CRMT = fielddata(:,12)/100;
delta_hHR = zeros(length(fielddata),1);
p1 = fielddata(:,2)*1e5;
pR = fielddata(:,15)*1e5;
xR = fielddata(:,24)/100;
TA0 = fielddata(:,3)+273.15;
ToMT = fielddata(:,13)+273.15;
TMTi = fielddata(:,11)+273.15;
% Declaration of vectors
dBP = NaN(length(fielddata),1);
dG = NaN(length(fielddata),1);
hG = NaN(length(fielddata),1);
hL = NaN(length(fielddata),1);
hHR = NaN(length(fielddata),1);
hHRd = NaN(length(fielddata),1);
DmQ = NaN(length(fielddata),1);
hBP = NaN(length(fielddata),1);
hR = NaN(length(fielddata),1);
qR = NaN(length(fielddata),1);
TA1 = NaN(length(fielddata),1);
h1 = NaN(length(fielddata),1);
pMT = NaN(length(fielddata),1);
dMT = NaN(length(fielddata),1);
% Conversions
for it = 1:length(fielddata)
    dG(it,1) = CoolProp.PropsSI('D','P',pR(it,1),'Q',1,'CO2');
    hG(it,1) = CoolProp.PropsSI('H','P',pR(it,1),'Q',1,'CO2');
    hL(it,1) = CoolProp.PropsSI('H','P',pR(it,1),'Q',0,'CO2');
    try
        hBP(it,1) = CoolProp.PropsSI('H','P',p1(it,1),'T',TBP(it,1),'CO2');
        dBP(it,1) = CoolProp.PropsSI('D','P',p1(it,1),'T',TBP(it,1),'CO2');
        if hBP(it) > CoolProp.PropsSI('H','P',p1(it,1),'Q',0,'CO2') % Right side of saturation curve
            hBP(it,1) = CoolProp.PropsSI('H','P',p1(it,1),'Q',0,'CO2') - 1e3;
            dBP(it,1) = CoolProp.PropsSI('D','P',p1(it,1),'Q',0,'CO2');
        end
    catch % Under saturation curve
        hBP(it,1) = CoolProp.PropsSI('H','P',p1(it,1),'Q',0,'CO2') - 1e3; % Substraction for safety
        dBP(it,1) = CoolProp.PropsSI('D','P',p1(it,1),'Q',0,'CO2');
    end
    hHR(it,1) = CoolProp.PropsSI('H','P',p1(it,1),'T',TMT(it,1),'CO2'); % No IT used.
    hHRd(it,1) = CoolProp.PropsSI('H','P',p1(it,1),'T',TMTd(it,1),'CO2'); % No IT used. Delay.
    qR(it,1) = x2q(xR(it,1),pR(it,1));
    hR(it,1) = CoolProp.PropsSI('H','P',pR(it,1),'Q',qR(it,1),'CO2');
    pMT(it,1) = CoolProp.PropsSI('P','T',ToMT(it,1),'Q',1,'CO2');
    dMT(it,1) = CoolProp.PropsSI('D','P',pMT(it,1),'T',TMTi(it,1),'CO2');
    h1(it,1) = 1/3*hHRd(it,1) + 2/3*hBP(it,1);
    DmQ(it,1) = dMT(it,1)*Mxf*VMT*CRMT(it,1);
    if it < TA0offsetindex
        TA0(it,1) = TA0(it,1) + 5;
    end
end
Y = [p1 hBP pR hR h1];
U = [CRA CRV CRIT CRG DmQ dBP dG hG hL hHR pMT TA0];

save('fielddata','fielddata','U','Y')
