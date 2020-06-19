#!/usr/bin/env wolframscript
(* ::Package:: *)

CloseAll;


F[z_]=a*z^2+b*z+c; (* Equazione della parabola *)
delta=b^2-4*a*c;
(* SISTEMA D'EQUAZIONI PER TROVARE LA FUNZIONE DELLA NOSTRA PARABOLA *)
E1=-b/(2*a)==l/2         (* Equazione per cordinata x del vertice della parabloa *)                    
E2=-delta/(4*a)==h     (* Equazione per cordinata y del vertice della parabloa *)          
E3=F[0]==0;                            (* Parabola passante per punto A ''il nostro vincolo'' *)            
S=Solve[{ E1, E2, E3},{a, b,c}]
G[z]=Simplify[F[z]/.{S}]


D1=Simplify[D[G[z]/.{S},z]] (* Derivata prima della funzione della parabola *)
D2=D[D1,z]  (* Derivata seconda della funzione della parabola *)


(* PARAMETRIZZAZIONE DELLA CURVA *)


 TAN\[Theta][z_]=D1             (* porto tang\[Theta] in funzione di z *)
ds=Simplify[Sqrt[1+D1^2]]       (* porto ds in funzione di z *)
SIN\[Theta][z_] =D1/ds       (* porto sin\[Theta] in funzione di z *)
COS\[Theta][z_] =1/ds          (* porto cos\[Theta] in funzione di z *)


(* SISTEMA FO *)


E4=(VA-VC)* COS\[Theta][0]-HA *SIN\[Theta][0]==0 ;          (* equazione alla traslazione orizzontale *)  
E5=l*p-HA *COS\[Theta][0]-(VA +VC)*SIN\[Theta][0]==0 ;       (* equazione alla traslazione verticale *)  
E6=l*VC *SIN\[Theta][0]-(l^2)*p/2==0 ;                      (* equazione alla rotazione in A *)  
SO=Solve[{E4, E5, E6},{VA,VC,HA}];\[NonBreakingSpace]\[NonBreakingSpace]


{HA,VA,VC}/. SO;


VAO={VA}/. SO
VCO={VC}/. SO 
HAO={HA}/. SO


M0[z_]=Simplify[{{{{-((p z^2)/2)+VAO*SIN\[Theta][0]*(z)-VAO*COS\[Theta][0]*G[z]}}}}]  (* Momento in F0 *)


(* SISTEMA F1 *) 


E7=VA*COS\[Theta] [0]-HA*SIN\[Theta][0]-VC*COS\[Theta][0] +1*SIN\[Theta][0]==0;       (* equazione alla traslazione orizzontale *)
E8=-VA*SIN\[Theta][0]-VC*SIN\[Theta][0]-HA*COS\[Theta][0]-1*COS\[Theta][0]==0;         (* equazione alla traslazione verticale *)
E9=VC*SIN\[Theta][0]*l+1*COS\[Theta][0]*l==0;                                 (* equazione alla rotazione in A *)   


S1=Solve[{ E7, E8, E9},{VA,VC,HA}];  


{HA,VA,VC}/. S1;


VA1={VA}/. S1
VC1={VC}/. S1 
HA1={HA}/. S1


M1[z_]=Simplify[{{{{+VAO*SIN\[Theta][0]*(z)-VAO*COS\[Theta][0]*G[z]+HA1*COS\[Theta][0]*(z)+HA1*SIN\[Theta][0]*G[z]}}}}]  (* Momento in F1 *) 


(* MULLER-BRESLAU *)
\[ScriptCapitalL]v10=Integrate[M1[z_]*M0[z_]/(EY*J)*ds,{z,0,l}]
\[ScriptCapitalL]v11 =Integrate[M1[z_]^2/(EY*J)*ds,{z,0,l}]
 \[Eta]10 =\[ScriptCapitalL]v10 ; \[Eta]11 =\[ScriptCapitalL]v11 ; 
\[Eta]1=\[Eta]10 +X1*\[Eta]11  ==0;
S2=Solve[{ \[Eta]1},{X1}]


X1=0  (* X1=0 perch\[EGrave]  \[Eta]10=0 ed  \[Eta]1=0 quindi anche senza determinare l'integrale de  \[Eta]11 sappiamo che X1 dev'essere nullo *)


(* Quindi le nostre reazioni effettive sono quelle nel sistema zero *)


VAE=VAO


VCE=VAO


HAE=HAO
