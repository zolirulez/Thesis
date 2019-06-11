tab = readtable('Gasloop_20190606.Csv');
tab.Properties.VariableNames'
figure(1)
plot(tab{:,2:end})
fielddata = tab{:,1.5e4:end};
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
CRA = fielddata(:,7)/100;
BP = zeros(length(fielddata),1);
CRV = fielddata(:,8)/100;
CRIT = fielddata(:,21)/100;
delta_hHR = zeros(length(fielddata),1);
p1 = fielddata(:,2)/100;
for it = 1:length(fielddata)
dBP = 
dG = 
hG =
hL =
hHR = 
DmQ =
hBP = 
pR = 
hR =
TA1 =
end
 
TA0 = fielddata(:,3)/100;
 
U = [CRA BP CRV CRIT delta_hHR dBP dG hG hL TA0 hHR DmQ];

 
% UFunction
%        CRA
%         BP
%        CRV
%       CRIT
%  delta_hHR
%        dBP
%         dG
%         hG
%         hL
%        TA0
%        hHR
%        DmQ
% YFunction
%   p1
%   hBP
%   pR
%   hR
%   TA1

