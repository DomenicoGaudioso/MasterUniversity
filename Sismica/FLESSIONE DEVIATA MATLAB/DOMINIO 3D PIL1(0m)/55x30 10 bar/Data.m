% DATI DI INPUT:

% - SEZIONE

% nome sezione
sezname = '55x30 10 bar';

% base sezione [mm]
bsez = 300;

% altezza sezione [mm]
hsez = 550;

% matrice barre d'armatura coordinate x-y e area [x1 y1 A1 ; .. .. .. ; x2 y2 a2] [mm][mm][mm^2]
XYAArm = [40 40 354 ; 150 40 254 ; 260 40 354 ; 40 183 254;260 183 254 ; 40 366 254 ; 260 366 254 ;40 510 254 ; 150 510 354 ; 260 510 354];


% - MATERIALI

% CLS

% valore resistenza a compressione [MPa]
fcd = 17;

% valore di deformazione limite elastico
e2 = 0.2;

% valore di deformazione ultima
eu = 0.35;

% ACCIAIO

% valore resistenza a trazione-compressione  [MPa]
fyd = 391.3;

% valore di deformazione limite elastico
ey = 0.1957;

% valore di deformazione ultima
esu = 6.75;


% - PARAMETRI PER CALCOLO DOMINIO

% vettore angoli di rotazione sezione [° sessa-decimali, da 0° a 90°]
angle = [0 15 30 45 60 75 90];

% valore di avanzamento (0.01 - 0.1)
step = 0.05;

% numero piani dominio 3D
slice = 40;
