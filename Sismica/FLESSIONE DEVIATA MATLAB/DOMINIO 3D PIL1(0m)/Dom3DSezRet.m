close all
clear all
clc

% Caricamento dati
Data;

% Creazione cartella di salvataggio
folder = [sezname];
mkdir(folder);

% copia file di input 'Data'
copyfile('Data.m',['./' sezname '/'])

% plottaggio sezione
SezRetPlot(bsez,hsez,XYAArm,sezname)


% CURVE MATERIALI
% Calcolo punti curva tensione-deformazione parabola rettangolo cls
[ecls,scls] = ParRetcls(fcd,e2,eu,sezname);
% Calcolo punti curva tensione-deformazione triangolo rettangolo acciaio
[eacc,sacc] = TriRetacc(fyd,ey,esu,sezname);


% 1° CICLO

% Dimensioni sezione
b = bsez;
h = hsez;

% Calcolo baricentro sezione
Xg = b/2;
Yg = h/2;

% Calcolo coordinate aree elementari riferite alla condizione iniziale
[XSezrifRot,YSezrifRot,Lrif,ANGrifRot,XSez,YSez,L,ANGRot] = SezRetRot(b,h,0);


for j=1:length(angle)

% Calcolo coordinate sezione ruotata
[XSezrifRot,YSezrifRot,Lrif,ANGrifRot,XSezRot,YSezRot,L,ANGRot] = SezRetRot(b,h,angle(j));
% Calcolo coordinate barre di armatura ruotate
[XYAArmRot,LArm,ANGArmRot] = ArmRot(XYAArm,angle(j));


% DEFINIZIONE CAMPI DI ROTTURA E CALCOLO COPPIE N-Mrd

% CAMPO 0 (trazione pura)  ||  ec=-esu  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], (corvatura) chi, N [kN], Mrd[kNm]
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,-esu,-esu);
[N01,M0x1,M0y1] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
% definizione limite campo di rottura
Nlim1 = N01;

% CAMPO 1 (tensoflessione)  ||  -esu<ec<=0  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], (corvatura) chi, N [kN], Mrd[kNm]
ec1 = [-esu+step:step:0];
for i=1:length(ec1)
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec1(i),-esu);
[N11(i),M1x1(i),M1y1(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim1 = [Nlim1 N11(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec1)
x11(i) = interp1([ec1(i),-esu],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 2 (trazione e compressione)  ||  0<ec<=eu  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
ec2 = [0+step:step:eu];
for i=1:length(ec2);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec2(i),-esu);
[N21(i),M2x1(i),M2y1(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim1 = [Nlim1 N21(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec2)
x21(i) = interp1([ec2(i),-esu],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 3 (trazione e compressione)  ||  ec=eu  -esu<es<=-ey
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
es3 = [-esu+step:step:-ey];
for i=1:length(es3);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es3(i));
[N31(i),M3x1(i),M3y1(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim1 = [Nlim1 N31(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es3)
x31(i) = interp1([eu,es3(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 4 (trazione e compressione)  ||  ec=eu  -ey<es<=0
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
es4 = [-ey+step:step:0];
for i=1:length(es4);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es4(i));
[N41(i),M4x1(i),M4y1(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim1 = [Nlim1 N41(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es4)
x41(i) = interp1([eu,es4(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 5 (pressoflessione)  ||  ec=eu  0<es<=esc
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
[esc] = escomp (YSezrifRot,XYAArmRot,eu);
es5 = linspace(0,esc,10);
for i=1:length(es5);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es5(i));
[N51(i),M5x1(i),M5y1(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim1 = [Nlim1 N51(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es5)
x51(i) = interp1([eu,es5(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 6 (pressoflessione)  ||  e2<=ec<eu  esc<es<=e2
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
ec6 = linspace(eu,e2,10);
es6 = linspace(esc,e2,10);
for i=1:length(ec6)-1;
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec6(i),es6(i));
[N61(i),M6x1(i),M6y1(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec6)-1
x61(i) = interp1([ec6(i),es6(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 7 (compressione pura)  ||  ec=e2  es=e2
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,e2,e2);
[N71,M7x1,M7y1] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
% definizione limite campo di rottura
Nlim1 = [Nlim1 N71];

% Salvataggio coppie N-Mrd
N1(:,j) = [N01 N11 N21 N31 N41 N51 N61 N71];
Mx1(:,j) = [M0x1 M1x1 M2x1 M3x1 M4x1 M5x1 M6x1 M7x1];
My1(:,j) = [M0y1 M1y1 M2y1 M3y1 M4y1 M5y1 M6y1 M7y1];
x1 (:,j) = [x11(1) x11 x21 x31 x41 x51 x61 x61(9)];
Nl1(:,j) = Nlim1;

end

% Definizione matrice armature per ciclo successivo
[ra ca] = size(XYAArmRot);
for i=1:ra
XYAArm(i,1) = h+XYAArmRot(i,1);
XYAArm(i,2) = XYAArmRot(i,2);
XYAArm(i,3) = XYAArmRot(i,3);
end




% 2° CICLO

% Dimensioni sezione
b = hsez;
h = bsez;

% Calcolo baricentro sezione
Xg = b/2;
Yg = h/2;

% Calcolo coordinate aree elementari riferite alla condizione iniziale
[XSezrifRot,YSezrifRot,Lrif,ANGrifRot,XSez,YSez,L,ANGRot] = SezRetRot(b,h,0);


for j=1:length(angle)

% Calcolo coordinate sezione ruotata
[XSezrifRot,YSezrifRot,Lrif,ANGrifRot,XSezRot,YSezRot,L,ANGRot] = SezRetRot(b,h,angle(j));
% Calcolo coordinate barre di armatura ruotate
[XYAArmRot,LArm,ANGArmRot] = ArmRot(XYAArm,angle(j));


% DEFINIZIONE CAMPI DI ROTTURA E CALCOLO COPPIE N-Mrd

% CAMPO 0 (trazione pura)  ||  ec=-esu  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], (corvatura) chi, N [kN], Mrd[kNm]
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,-esu,-esu);
[N02,M0y2,M0x2] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
% definizione limte campo di rottura
Nlim2 = N02;

% CAMPO 1 (tensoflessione)  ||  -esu<ec<=0  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], (corvatura) chi, N [kN], Mrd[kNm]
ec1 = [-esu+step:step:0];
for i=1:length(ec1)
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec1(i),-esu);
[N12(i),M1y2(i),M1x2(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim2 = [Nlim2 N12(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec1)
x12(i) = interp1([ec1(i),-esu],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 2 (trazione e compressione)  ||  0<ec<=eu  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
ec2 = [0+step:step:eu];
for i=1:length(ec2);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec2(i),-esu);
[N22(i),M2y2(i),M2x2(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim2 = [Nlim2 N22(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec2)
x22(i) = interp1([ec2(i),-esu],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 3 (trazione e compressione)  ||  ec=eu  -esu<es<=-ey
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
es3 = [-esu+step:step:-ey];
for i=1:length(es3);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es3(i));
[N32(i),M3y2(i),M3x2(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim2 = [Nlim2 N32(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es3)
x32(i) = interp1([eu,es3(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 4 (trazione e compressione)  ||  ec=eu  -ey<es<=0
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
es4 = [-ey+step:step:0];
for i=1:length(es4);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es4(i));
[N42(i),M4y2(i),M4x2(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim2 = [Nlim2 N42(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es4)
x42(i) = interp1([eu,es4(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 5 (pressoflessione)  ||  ec=eu  0<es<=esc
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
[esc] = escomp (YSezrifRot,XYAArmRot,eu);
es5 = linspace(0,esc,10);
for i=1:length(es5);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es5(i));
[N52(i),M5y2(i),M5x2(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim2 = [Nlim2 N52(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es5)
x52(i) = interp1([eu,es5(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 6 (pressoflessione)  ||  e2<=ec<eu  esc<es<=e2
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
ec6 = linspace(eu,e2,10);
es6 = linspace(esc,e2,10);
for i=1:length(ec6)-1;
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec6(i),es6(i));
[N62(i),M6y2(i),M6x2(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec6)-1
x62(i) = interp1([ec6(i),es6(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 7 (compressione pura)  ||  ec=e2  es=e2
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,e2,e2);
[N72,M7y2,M7x2] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
% definizione limite campo di rottura
Nlim2 = [Nlim2 N72];

% Salvataggio coppie N-Mrd
N2(:,j) = [N02 N12 N22 N32 N42 N52 N62 N72];
Mx2(:,j) = -1*[M0x2 M1x2 M2x2 M3x2 M4x2 M5x2 M6x2 M7x2];
My2(:,j) = [M0y2 M1y2 M2y2 M3y2 M4y2 M5y2 M6y2 M7y2];
x2 (:,j) = [x12(1) x12 x22 x32 x42 x52 x62 x62(9)];
Nl2 (:,j) = Nlim2;

end

% Definizione matrice armature per ciclo successivo
for i=1:ra
XYAArm(i,1) = h+XYAArmRot(i,1);
XYAArm(i,2) = XYAArmRot(i,2);
XYAArm(i,3) = XYAArmRot(i,3);
end




% 3° CICLO

% Dimensioni sezione
b = bsez;
h = hsez;

% Calcolo baricentro sezione
Xg = b/2;
Yg = h/2;

% Calcolo coordinate aree elementari riferite alla condizione iniziale
[XSezrifRot,YSezrifRot,Lrif,ANGrifRot,XSez,YSez,L,ANGRot] = SezRetRot(b,h,0);


for j=1:length(angle)

% Calcolo coordinate sezione ruotata
[XSezrifRot,YSezrifRot,Lrif,ANGrifRot,XSezRot,YSezRot,L,ANGRot] = SezRetRot(b,h,angle(j));
% Calcolo coordinate barre di armatura ruotate
[XYAArmRot,LArm,ANGArmRot] = ArmRot(XYAArm,angle(j));


% DEFINIZIONE CAMPI DI ROTTURA E CALCOLO COPPIE N-Mrd

% CAMPO 0 (trazione pura)  ||  ec=-esu  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], (corvatura) chi, N [kN], Mrd[kNm]
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,-esu,-esu);
[N03,M0x3,M0y3] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
% definizione limite campo di rottura
Nlim3 = N03;

% CAMPO 1 (tensoflessione)  ||  -esu<ec<=0  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], (corvatura) chi, N [kN], Mrd[kNm]
ec1 = [-esu+step:step:0];
for i=1:length(ec1)
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec1(i),-esu);
[N13(i),M1x3(i),M1y3(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim3 = [Nlim3 N13(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec1)
x13(i) = interp1([ec1(i),-esu],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 2 (trazione e compressione)  ||  0<ec<=eu  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
ec2 = [0+step:step:eu];
for i=1:length(ec2);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec2(i),-esu);
[N23(i),M2x3(i),M2y3(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim3 = [Nlim3 N23(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec2)
x23(i) = interp1([ec2(i),-esu],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 3 (trazione e compressione)  ||  ec=eu  -esu<es<=-ey
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
es3 = [-esu+step:step:-ey];
for i=1:length(es3);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es3(i));
[N33(i),M3x3(i),M3y3(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim3 = [Nlim3 N33(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es3)
x33(i) = interp1([eu,es3(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 4 (trazione e compressione)  ||  ec=eu  -ey<es<=0
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
es4 = [-ey+step:step:0];
for i=1:length(es4);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es4(i));
[N43(i),M4x3(i),M4y3(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim3 = [Nlim3 N43(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es4)
x43(i) = interp1([eu,es4(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 5 (pressoflessione)  ||  ec=eu  0<es<=esc
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
[esc] = escomp (YSezrifRot,XYAArmRot,eu);
es5 = linspace(0,esc,10);
for i=1:length(es5);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es5(i));
[N53(i),M5x3(i),M5y3(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim3 = [Nlim3 N53(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es5)
x53(i) = interp1([eu,es5(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 6 (pressoflessione)  ||  e2<=ec<eu  esc<es<=e2
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
ec6 = linspace(eu,e2,10);
es6 = linspace(esc,e2,10);
for i=1:length(ec6)-1;
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec6(i),es6(i));
[N63(i),M6x3(i),M6y3(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec6)-1
x63(i) = interp1([ec6(i),es6(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 7 (compressione pura)  ||  ec=e2  es=e2
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,e2,e2);
[N73,M7x3,M7y3] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
% definizione limite campo di rottura
Nlim3 = [Nlim3 N73];

% Salvataggio coppie N-Mrd
N3(:,j) = [N03 N13 N23 N33 N43 N53 N63 N73];
Mx3(:,j) = -1*[M0x3 M1x3 M2x3 M3x3 M4x3 M5x3 M6x3 M7x3];
My3(:,j) = -1*[ M0y3 M1y3 M2y3 M3y3 M4y3 M5y3 M6y3 M7y3];
x3 (:,j) = [x13(1) x13 x23 x33 x43 x53 x63 x63(9)];
Nl3 (:,j) = Nlim3;

end

% Definizione matrice armature per ciclo successivo
for i=1:ra
XYAArm(i,1) = h+XYAArmRot(i,1);
XYAArm(i,2) = XYAArmRot(i,2);
XYAArm(i,3) = XYAArmRot(i,3);
end




% 4° CICLO

% Dimensioni sezione
b = hsez;
h = bsez;

% Calcolo baricentro sezione
Xg = b/2;
Yg = h/2;

% Calcolo coordinate aree elementari riferite alla condizione iniziale
[XSezrifRot,YSezrifRot,Lrif,ANGrifRot,XSez,YSez,L,ANGRot] = SezRetRot(b,h,0);


for j=1:length(angle)

% Calcolo coordinate sezione ruotata
[XSezrifRot,YSezrifRot,Lrif,ANGrifRot,XSezRot,YSezRot,L,ANGRot] = SezRetRot(b,h,angle(j));
% Calcolo coordinate barre di armatura ruotate
[XYAArmRot,LArm,ANGArmRot] = ArmRot(XYAArm,angle(j));


% DEFINIZIONE CAMPI DI ROTTURA E CALCOLO COPPIE N-Mrd

% CAMPO 0 (trazione pura)  ||  ec=-esu  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], (corvatura) chi, N [kN], Mrd[kNm]
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,-esu,-esu);
[N04,M0y4,M0x4] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
% definizione limite campo di rottura
Nlim4 = N04;

% CAMPO 1 (tensoflessione)  ||  -esu<ec<=0  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], (corvatura) chi, N [kN], Mrd[kNm]
ec1 = [-esu+step:step:0];
for i=1:length(ec1)
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec1(i),-esu);
[N14(i),M1y4(i),M1x4(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim4 = [Nlim4 N14(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec1)
x14(i) = interp1([ec1(i),-esu],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 2 (trazione e compressione)  ||  0<ec<=eu  es=-esu
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
ec2 = [0+step:step:eu];
for i=1:length(ec2);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec2(i),-esu);
[N24(i),M2y4(i),M2x4(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim4 = [Nlim4 N24(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec2)
x24(i) = interp1([ec2(i),-esu],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 3 (trazione e compressione)  ||  ec=eu  -esu<es<=-ey
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
es3 = [-esu+step:step:-ey];
for i=1:length(es3);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es3(i));
[N34(i),M3y4(i),M3x4(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim4 = [Nlim4 N34(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es3)
x34(i) = interp1([eu,es3(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 4 (trazione e compressione)  ||  ec=eu  -ey<es<=0
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
es4 = [-ey+step:step:0];
for i=1:length(es4);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es4(i));
[N44(i),M4y4(i),M4x4(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim4 = [Nlim4 N44(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es4)
x44(i) = interp1([eu,es4(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 5 (pressoflessione)  ||  ec=eu  0<es<=esc
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
[esc] = escomp (YSezrifRot,XYAArmRot,eu);
es5 = linspace(0,esc,10);
for i=1:length(es5);
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,eu,es5(i));
[N54(i),M5y4(i),M5x4(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% definizione limite campo di rottura
Nlim4 = [Nlim4 N54(i)];
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(es5)
x54(i) = interp1([eu,es5(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 6 (pressoflessione)  ||  e2<=ec<eu  esc<es<=e2
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
ec6 = linspace(eu,e2,10);
es6 = linspace(esc,e2,10);
for i=1:length(ec6)-1;
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,ec6(i),es6(i));
[N64(i),M6y4(i),M6x4(i)] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
end
% calcolo posizione x asse neutro dalla fibra superiore
for i=1:length(ec6)-1
x64(i) = interp1([ec6(i),es6(i)],[max(max(YSezrifRot)), min(XYAArmRot(:,2))],0,'linear','extrap');
end

% CAMPO 7 (compressione pura)  ||  ec=e2  es=e2
% Calcolo deformazioni sezione cls e barre d'armatura,(posizione asse neutro) x [mm], N [kN], Mrd[kNm]
[DefSez,DefArm] = SezRetDef(YSezRot,YSezrifRot,XYAArmRot,e2,e2);
[N74,M7y4,M7x4] = Campon (Xg,Yg,XSez,YSez,XYAArm,DefSez,DefArm,ecls,scls,eacc,sacc);
% definizione limite campo di rottura
Nlim4 = [Nlim4 N74];

% Salvataggio coppie N-Mrd
N4(:,j) = [N04 N14 N24 N34 N44 N54 N64 N74];
Mx4(:,j) = [M0x4 M1x4 M2x4 M3x4 M4x4 M5x4 M6x4 M7x4];
My4(:,j) = -1*[M0y4 M1y4 M2y4 M3y4 M4y4 M5y4 M6y4 M7y4];
x4 (:,j) = [x14(1) x14 x24 x34 x44 x54 x64 x64(9)];
Nl4 (:,j) = Nlim4;

end



% assemblaggio matrici dominio di resistenza
N = [N1 N2 N3 N4];
Mx = [Mx1 Mx2 Mx3 Mx4];
My = [My1 My2 My3 My4];
x = [x1 x2 x3 x4];
Nl = [Nl1 Nl2 Nl3 Nl4];

% salvataggio matrici dominio di resistenza
folder = ['./' sezname '/DATI DOMINIO/COMPLETI'];
mkdir(folder)
filename = [folder '/N.txt'];
dlmwrite(filename,N,'delimiter','\t','precision',15)
filename = [folder '/Mx.txt'];
dlmwrite(filename,Mx,'delimiter','\t','precision',15)
filename = [folder '/My.txt'];
dlmwrite(filename,My,'delimiter','\t','precision',15)
filename = [folder '/x.txt'];
dlmwrite(filename,x,'delimiter','\t','precision',15)
% angolo asse neutro per ogni colonna di N Mx My (si parte dalla posizione orizzontale ruotando in senso orario)
ANGn1 = [0:90/(length(angle)-1):90];
ANGn2 = [90:90/(length(angle)-1):180];
ANGn3 = [180:90/(length(angle)-1):270];
ANGn4 = [270:90/(length(angle)-1):360];
ANGn = [ANGn1 ANGn2 ANGn3 ANGn4];
filename = [folder '/ANGn.txt'];
dlmwrite(filename,ANGn,'delimiter','\t','precision',4)

% Plottaggio grafico campi di rottura
figure(250)
hold on
grid on
title(['LIMITI CAMPI DI ROTTURA "' sezname '"'])
xlabel('Angolo inclinazione asse neutro')
ylabel('Sforzo normale  [kN]')
xlim([0 360])
for i=1:7
plot(ANGn,Nl(i,:),'-o')
end
filename = ['./' sezname '/SEZIONE' sezname ' LIMITI CAMPI'];
saveas(gcf,filename,'png');

% Plottaggio dominio di resistenza
figure(300)
% grafico sollecitazioni
NMxMy = dlmread('NMxMy.txt');
for o=1:(length(NMxMy)/2)
pil1(:,:,o)=[NMxMy(2*o-1,:);NMxMy(2*o,:)];
scatter3(pil1(:,2,o),pil1(:,3,o),pil1(:,1,o),'filled')%'MarkerFaceColor','r'
view(80,20)
hold on
end

grid on
title(['DOMINIO DI RESISTENZA "' sezname '"'])
xlabel('Momento flettente Mx [kNm]')
ylabel('Momento flettente My [kNm]')
zlabel('Sforzo normale  [kN]')
plot3(Mx,My,N,'b')




% Plottaggio isolinee dominio di resistenza
% Calcolo matrici di interpolazione (la matrice delle azioni assiali viene ripulita dai valori ripetuti per consentire l'interpolazione)
[r c] = size(N);
for i=1:c
[Nclean(:,i)] = cleancurve (N(:,i));
end
% vettore di campionamento
Nval = linspace(max(min(Nclean)),min(max(Nclean)),slice);
% interpolazione
[r c] = size(Nclean);
for i=1:c
Mxval(:,i) = interp1(Nclean(:,i),Mx(:,i),Nval,'linear');
Myval(:,i) = interp1(Nclean(:,i),My(:,i),Nval,'linear');
xval (:,i) = interp1(Nclean(:,i),x(:,i),Nval,'linear');
end
% creazione vettore sforzo assiale costante (isolinea) e plottaggio
for p=1:slice
Nxy = Nval(p) * ones(1,c);
plot3(Mxval(p,:),Myval(p,:),Nxy,'r')
end
% salvataggio
filename = ['./' sezname '/SEZIONE' sezname ' DOMINIO 3D linee'];
saveas(gcf,filename,'png');


% Plottaggio superficie dominio di resistenza
figure(400)
hold on
grid on
title(['DOMINIO DI RESISTENZA "' sezname '"'])
xlabel('Momento flettente Mx [kNm]')
ylabel('Momento flettente My [kNm]')
zlabel('Sforzo normale  [kN]')
% creazione matrice degli sforzi assiali (tutte le colonne sono uguali)
for k=1:4*length(angle)
Nvalmatrix(:,k) = Nval;
end
surf(Mxval,Myval,Nvalmatrix)
view(40,20)
colorbar
filename = ['./' sezname '/SEZIONE' sezname ' DOMINIO 3D superficie'];
saveas(gcf,filename,'png');


% salvataggio variabili dominio 3D
filename = ['./' sezname '/DATI DOMINIO/N 3D.txt'];
dlmwrite(filename,Nvalmatrix,'delimiter','\t','precision',16)
filename = ['./' sezname '/DATI DOMINIO/Mx 3D.txt'];
dlmwrite(filename,Mxval,'delimiter','\t','precision',16)
filename = ['./' sezname '/DATI DOMINIO/My 3D.txt'];
dlmwrite(filename,Myval,'delimiter','\t','precision',16)
filename = ['./' sezname '/DATI DOMINIO/x 3D.txt'];
dlmwrite(filename,xval,'delimiter','\t','precision',16)


% verifica sezioni
folder = ['./' sezname '/VERIFICHE'];
mkdir(folder)
NMxMy = dlmread('NMxMy.txt');
[VER] = versezret (folder,Nval,Mxval,Myval,NMxMy);
% salvataggio dati verifiche
filename = [folder '/VERIFICHE.txt'];
dlmwrite(filename,VER,'delimiter','\t','precision',16)