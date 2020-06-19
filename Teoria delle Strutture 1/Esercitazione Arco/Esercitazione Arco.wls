#!/usr/bin/env wolframscript
(* ::Package:: *)

(* PUNTO 1 *)
f=L/2;
F[z_]=a*z^2+b*z+c;                 (* Equazione della parabola *)
delta=b^2-4*a*c;       (* SISTEMA D'EQUAZIONI PER TROVARE LA FUNZIONE DELLA PARABOLA DATA *)
E1=-b/(2*a)==L/2;         (* Equazione per cordinata x del vertice della parabloa *)                    
E2=-delta/(4*a)==f;       (* Equazione per cordinata y del vertice della parabloa *)          
E3=F[0] ==0;               (* Parabola passante per punto A ''il nostro vincolo'' *)            
S=Solve[{ E1, E2, E3},{a, b,c}]
G[z_]=Simplify[F[z]/.{S}]


D1=Simplify[D[G[z]/.{S},z]]  (* Derivata prima della funzione della parabola *)
D2=D[D1,z]                   (* Derivata seconda della funzione della parabola *)


(* PARAMETRIZZAZIONE DELLA CURVA *)
TAN\[Theta][z_]=D1                  (* porto tang\[Theta] in funzione di z *)
ds=Simplify[Sqrt[1+D1^2]]     (* porto ds in funzione di z *)
SIN\[Theta][z_] =D1/ds              (* porto sin\[Theta] in funzione di z *)
COS\[Theta][z_] =1/ds               (* porto cos\[Theta] in funzione di z *)


(* SISTEMA FO *)
(* equazione alla traslazione orizzontale *)  
E4=(VA-VC)* COS\[Theta][0]-HA *SIN\[Theta][0]==0 ;     
(* equazione alla traslazione verticale *)   
E5=L *p-HA *COS\[Theta][0]-(VA +VC)*SIN\[Theta][0]==0 ; 
(* equazione alla rotazione in A *)   
E6=L*VC *SIN\[Theta][0]-(L^2)*p/2==0 ;                
SO=Solve[{ E4, E5, E6},{VA,VC,HA}];\[NonBreakingSpace]
HAO={HA}/. SO 
VAO={VA}/. SO 
VCO={VC}/. SO \[NonBreakingSpace]


(* Momento in F0 *)
M0[z_]=Simplify[-((p z^2)/2)+VAO*SIN\[Theta][0]*(z)-VAO *COS\[Theta][0]*G[z] ]  


(* SISTEMA F1 *) 
(* equazione alla traslazione orizzontale *)
E7=VA*COS\[Theta] [0]-HA*SIN\[Theta][0]-VC*COS\[Theta][0]+1*SIN\[Theta][0]==0;   
(* equazione alla traslazione verticale *)      
E8=-VA*SIN\[Theta][0]-VC*SIN\[Theta][0]-HA*COS\[Theta][0]-1*COS\[Theta][0]==0;   
(* equazione alla rotazione in A *)       
E9=VC*SIN\[Theta][0]*L+1*COS\[Theta] [0]*L==0;                                 
S1=Solve[{ E7, E8, E9},{VA,VC,HA}]; 
HA1={HA}/. S1 
VA1={VA}/. S1 
VC1={VC}/. S1   


(* Momento in F1 *)
M1[z_]=Simplify[VA1*SIN\[Theta][0]*(z)-VA1*COS\[Theta][0]*G[z]+1*COS\[Theta][0]*(z)+1*SIN\[Theta][0]*(G[z])] 


(* MULLER-BRESLAU *)
\[ScriptCapitalL]v10=Simplify[ Integrate[M1[z]*M0[z]/(EY*J)*ds,{z,0,L}]];
\[ScriptCapitalL]v11 =Integrate[M1[z]^2/(EY*J)*ds,{z,0,L}];
 \[Eta]10 =\[ScriptCapitalL]v10 ; \[Eta]11 =\[ScriptCapitalL]v11 ; 
\[Eta]1=\[Eta]10 +X1*\[Eta]11  ==0;
S2=Solve[{ \[Eta]1},{X1}]
XE=0


M1E=M0[z]+XE*M1[z] (* Momento esatto del PUNTO 1 *)


VAE=VAO+XE*VA1  (* VA esatto del PUNTO 1 *)


HAE=HAO+XE*HA1\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace](* HA esatto del PUNTO 1 *) 


(* CALCOLO DELLE CDS *)
(* Forza normale di compressione del PUNTO 1 *)  
N1E=Simplify[ -VAE*(COS\[Theta][0])/COS\[Theta][z]]  


(* PUNTO 2 *)
(* equazione alla traslazione orizzontale *)\[NonBreakingSpace]\[NonBreakingSpace]
E10=(VA-VC)* COS\[Theta][0]-HA *SIN\[Theta][0]+q*L ==0 ;\[NonBreakingSpace]
(* equazione alla traslazione verticale *)\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]
E11=-HA *COS\[Theta][0]-(VA +VC)*SIN\[Theta][0]==0 ;\[NonBreakingSpace]
(* equazione alla rotazione in A *)\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]             
E12=L*VC *SIN\[Theta][0]-q*L*f==0 ;\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]  
S3=Solve[{ E10,E11,E12},{VA,VC,HA}];
HA2O={HA}/. S3 
VA2O={VA}/. S3 
VC2O={VC}/. S3  


(* Momento in F0 del Punto 2 *)
M20[z_]=VA2O*SIN\[Theta][0]*z-VA2O*COS\[Theta][0]*G[z]+HA2O*COS\[Theta][0]*z+HA2O*SIN\[Theta][0]*G[z]+q*z*(f-G[z])\[NonBreakingSpace] 


(* SISTEMA F1 *)  
(* equazione alla traslazione orizzontale *)
E13=VA*COS\[Theta] [0]-HA*SIN\[Theta][0]-VC*COS\[Theta][0]+1*SIN\[Theta][0]==0;\[NonBreakingSpace]
(* equazione alla traslazione verticale *)\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]
E14=-VA*SIN\[Theta][0]-VC*SIN\[Theta][0]-HA*COS\[Theta][0]-1*COS\[Theta][0]==0;
(* equazione alla rotazione in A *)\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]
E15=VC*SIN\[Theta][0]*L+1*COS\[Theta] [0]*L==0;\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace] 
S4=Solve[{ E13,E14,E15},{VA,VC,HA}];
HA21={HA}/. S4 
VA21={VA}/. S4 
VC21={VC}/. S4\[NonBreakingSpace]\[NonBreakingSpace]\[NonBreakingSpace]


(* Momento in F1 del PUNTO 2 *)
M21[z_]=Simplify[VA21*SIN\[Theta][0]*z-VA21*COS\[Theta][0]*G[z]+1*COS\[Theta][0]*z+1*SIN\[Theta][0]*G[z]]\[NonBreakingSpace]\[NonBreakingSpace]


(* MULLER-BRESLAU *) 
\[ScriptCapitalL]v210=Simplify[ Integrate[M21[z]*M20[z]/(EY*J)*ds,{z,0,L}]]
\[ScriptCapitalL]v211 =Simplify[Integrate[M21[z]^2/(EY*J)*ds,{z,0,L}]]
\[Eta]210 =\[ScriptCapitalL]v210 ; \[Eta]211 =\[ScriptCapitalL]v211 ; 
\[Eta]21=\[Eta]210 +X21*\[Eta]211==0;
S5=Solve[{ \[Eta]21},{X21}]


X2E=Simplify[{X21}/. S5 ] 
M2E=Simplify[M20[z]+Simplify[X2E*M21[z]]] (* M[z] esatto del PUNTO 2*)
VA2E=Simplify[VA2O+Simplify[X2E*VA21]]     (* VA esatto del PUNTO 2*)
HA2E=Simplify[HA2O+Simplify[X2E*HA21]]     (* HA esatto del PUNTO 2 *)


(* FORZA NORMALE E TAGLIO PUNTO 2 *)
E16=N2*COS\[Theta][z]+T2*SIN\[Theta][z]+VA2E*COS\[Theta][0]-HA2E*SIN\[Theta][0]+q*z==0;
E17=T2*COS\[Theta][z]-N2*SIN\[Theta][z]-HA2E*COS\[Theta][0]-VA2E*SIN\[Theta][0]==0;
S6=Solve[{ E16,E17},{N2,T2}];
N2E=Simplify[{N2}/. S6]
T2E=Simplify[{T2}/. S6]


(* SOVRAPPOSIZIONE DEGLI EFFETTI PUNTO 1 + PUNTO 2 *)
MF=M2E
NF=Simplify[N1E+N2E]
TF=Simplify[T2E]


(* PUNTO 3 ECCENTRICITA *)
e=MF/NF


(* Derivata prima della funzione ''Curva delle Pressioni'' *)
TANe\[Theta][z_]=Simplify[D[e,z]]  


(* funzione dell' intradosso della parabola *)
H=L/20; (* spessore arco *)
(* Equazione per cordinata x del vertice della parabloa *)  
E18=-b/(2*a)==L/2;
(* Equazione per cordinata y del vertice della parabloa *)                                                      
E19=-delta/(4*a)==f-H/2;  
(* Parabola passante per punto A ''il nostro vincolo'' *)                                       
E20=F[(H/2)*SIN\[Theta][0]] ==-(H/2)*COS\[Theta][0];                        
S7=Solve[{ E18, E19, E20},{a, b,c}]
G2[z]=Simplify[F[z]/.{S7}]


(* funzione dell' estradosso della parabola *)
(* Equazione per cordinata x del vertice della parabloa *)   
E21=-b/(2*a)==L/2;   
(* Equazione per cordinata y del vertice della parabloa *)                                                   
E22=-delta/(4*a)==f+H/2;
(* Parabola passante per punto A ''il nostro vincolo'' *)                                        
E23=F[-(H/2)*SIN\[Theta][0]] ==(H/2)*COS\[Theta][0];                        
S8=Solve[{ E21, E22, E23},{a, b,c}]
G3[z]=Simplify[F[z]/.{S8}]


(* Derivata prima della funzione che descrive l'intradosso *)
TANintradosso\[Theta][z_]=Simplify[D[G2[z],z]]  
(* Derivata prima della funzione che descrive l'estradosso *)
TANestradosso\[Theta][z_]=Simplify[D[G3[z],z]]  


(*valore di q per cui la parabola \[EGrave] tangente all'estradosso o intradosso dell'arco*)
E24=TANintradosso\[Theta][z]==TANe\[Theta][z] 
S9=Solve[{E24},{q}];
q2=Simplify[{q}/.{S9}];

(* questo non basta, deve essere uguale anche fissando i punti dell'estradosso o intradosso*)
(* \[EGrave] IMPOSSIBILE che la curva sia tangente in tutti i punti dell'intradosso o estradosso *)
