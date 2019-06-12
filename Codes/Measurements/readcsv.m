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
VR = 0.155;
% Data
CRA = fielddata(:,7)/100;
TBP = fielddata(:,5)+273.15;
TMT = fielddata(:,10)+273.15;
BP = zeros(length(fielddata),1);
CRV = fielddata(:,8)/100;
CRIT = fielddata(:,21)/100;
delta_hHR = zeros(length(fielddata),1);
p1 = fielddata(:,2)*1e5;
pR = fielddata(:,15)*1e5;
xR = fielddata(:,24)/100;
TA0 = fielddata(:,3)+273.15;
dBP = NaN(length(fielddata),1);
dG = NaN(length(fielddata),1);
hG = NaN(length(fielddata),1);
hL = NaN(length(fielddata),1);
hHR = NaN(length(fielddata),1);
DmQ = NaN(length(fielddata),1);
hBP = NaN(length(fielddata),1);
hR = NaN(length(fielddata),1);
qR = NaN(length(fielddata),1);
TA1 = NaN(length(fielddata),1);
delay = 300;
for it = 1:length(fielddata)
    try
        dBP(it,1) = CoolProp.PropsSI('D','P',p1(it,1),'T',TBP(it,1),'CO2');
    catch
        dBP(it,1) = CoolProp.PropsSI('D','P',p1(it,1),'Q',0,'CO2');
    end
    dG(it,1) = CoolProp.PropsSI('D','P',pR(it,1),'Q',1,'CO2');
    hG(it,1) = CoolProp.PropsSI('H','P',pR(it,1),'Q',1,'CO2');
    hL(it,1) = CoolProp.PropsSI('H','P',pR(it,1),'Q',0,'CO2');
    try
        hBP(it,1) = CoolProp.PropsSI('H','P',p1(it,1),'T',TBP(it,1),'CO2');
    catch
        hBP(it,1) = CoolProp.PropsSI('H','P',p1(it,1),'Q',0,'CO2');
    end
    hHR(it,1) = CoolProp.PropsSI('H','P',p1(it,1),'T',TMT(it,1),'CO2'); % Not IT used.
    qR(it,1) = x2q(xR(it,1),pR(it,1));
    hR(it,1) = CoolProp.PropsSI('H','P',pR(it,1),'Q',qR(it,1),'CO2');
    DmQ(it,1) = 0;
    h1(it,1) = 1e5*(1.327157688246472+...
        0.003577339990579*(tab{it+start-1-delay,10}+273.15)+...
        0.003431954477368*(tab{it+start-1,3}+273.15)+...
        0.987387385770585*(tab{it+start-1-delay,19}/100)-...
        3.177645289705887*(tab{it+start-1-delay,7}/100));
end
U = [CRA BP CRV CRIT delta_hHR dBP dG hG hL TA0 hHR DmQ];
Y = [p1 hBP pR hR h1];
 
