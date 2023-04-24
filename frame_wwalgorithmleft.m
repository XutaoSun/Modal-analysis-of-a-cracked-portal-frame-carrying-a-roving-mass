function[WW]=frame_WWalgorithmLeft(w,L,nodes,conn,angle,sup1,sup2,J,M,E,G,k,I,Ip,A,Rho,Mu,Lcr11,Lcr22,Lcr33,Mcr11,Mcr22,Mcr33)
%Code to be stored in a file named 'WWalgorithmLeft.m'. 'Left' represents that the spring-mass system is roving on the left segment of the cracked beam
%'w' is the trial frequency in rad/s.
DFGstiff=zeros(3*nodes,3*nodes); %Global dynamic stiffness matrix of the cracked beam.
DFstiff=zeros(6,6); %dynamic stiffness matrix of each element, including the crack element.
Lcrackstiff=[1/Lcr11 0 0 -1/Lcr11 0 0; 0 1/Lcr22 0 0 -1/Lcr22 0; 0 0 1/Lcr33 0 0 -1/Lcr33; -1/Lcr11 0 0 1/Lcr11 0 0; 0 -1/Lcr22 0 0 1/Lcr22 0; 0 0 -1/Lcr33 0 0 1/Lcr33]; 
Mcrackstiff=[1/Mcr11 0 0 -1/Mcr11 0 0; 0 1/Mcr22 0 0 -1/Mcr22 0; 0 0 1/Mcr33 0 0 -1/Mcr33; -1/Mcr11 0 0 1/Mcr11 0 0; 0 -1/Mcr22 0 0 1/Mcr22 0; 0 0 -1/Mcr33 0 0 1/Mcr33];
Jm=0; ng=0; %Initial values for the terms of the Wittrick-Williams algorithm, where 'Jm' is the number of fixed end natural frequencies below the trial frequency 'w', 'ng' is the number of negative leading diagonal elements of the upper triangular matrix formed from the global stiffness matrix after performing Gaussian elimination without row interchange.
for i=1
	Alpha2(i)=Rho*A*w^2*L(i)^2/(E*A);
	Beta2(i)=Rho*Ip*Mu^2*w^2/(E*A);
	Gamma(i)=(Alpha2(i)/(1-Beta2(i)))^0.5;
	bb(i)=Rho*A*w^2*L(i)^4/(E*I);
	rr(i)=I/(A*L(i)^2);
	ss(i)=E*I/(k*A*G*L(i)^2);
	a1(i)=E*A/L(i)*Gamma(i)*(1-Beta2(i))*cot(Gamma(i));
	a2(i)=-E*A/L(i)*Gamma(i)*(1-Beta2(i))*csc(Gamma(i));
	Phi(i)=(bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	S(i)=sin(Phi(i));
	C(i)=cos(Phi(i));
	W1(i)=E*I/L(i);
	W2(i)=E*I/L(i)^2;
	W3(i)=E*I/L(i)^3;
	Z(i)=Phi(i)-bb(i)*ss(i)/Phi(i);
	if bb(i)*rr(i)*ss(i)<1
		j=1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sinh(Lambda(i));
		C_(i)=cosh(Lambda(i));	
	end
	if bb(i)*rr(i)*ss(i)>1
		j=-1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sin(Lambda(i));
		C_(i)=cos(Lambda(i));	
	end
	Eta(i)=Z(i)/(j*Lambda(i)+bb(i)*ss(i)/Lambda(i));
	Tau(i)=(Lambda(i)+Eta(i)*Phi(i))/(2*Eta(i)*(1-C(i)*C_(i))+(1-j*Eta(i)^2)*S(i)*S_(i));
	d1(i)=W3(i)*bb(i)*Tau(i)*(C(i)*S_(i)+Eta(i)*S(i)*C_(i))/(Lambda(i)*Phi(i));
	d2(i)=W2(i)*Z(i)*Tau(i)*((Phi(i)+j*Eta(i)*Lambda(i))*S(i)*S_(i)-(Lambda(i)-Eta(i)*Phi(i))*(1-C(i)*C_(i)))/(Lambda(i)+Eta(i)*Phi(i));
	d3(i)=W1(i)*Tau(i)*(S(i)*C_(i)-j*Eta(i)*C(i)*S_(i));
	d4(i)=-W3(i)*bb(i)*Tau(i)*(S_(i)+Eta(i)*S(i))/(Lambda(i)*Phi(i));
	d5(i)=W2(i)*Z(i)*Tau(i)*(C_(i)-C(i));
	d6(i)=W1(i)*Tau(i)*(j*Eta(i)*S_(i)-S(i));
	Jm=Jm+int16(fix(Gamma(i)/pi));
	DFstiff=[a1(i) 0 0 a2(i) 0 0; 0 d1(i) d2(i) 0 d4(i) d5(i); 0 d2(i) d3(i) 0 -d5(i) d6(i); a2(i) 0 0 a1(i) 0 0; 0 d4(i) -d5(i) 0 d1(i) -d2(i); 0 d5(i) d6(i) 0 -d2(i) d3(i)];
	t=[cos(angle(i)) sin(angle(i)) 0; -sin(angle(i)) cos(angle(i)) 0; 0 0 1];
	T=[t zeros(3,3); zeros(3,3) t];
	DFstiff=T'*DFstiff*T;
	n1=conn(i,1);
	n2=conn(i,2);
	DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)=DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)+DFstiff(1:3,1:3);
	DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)=DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)+DFstiff(1:3,4:6);
	DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)=DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)+DFstiff(4:6,1:3);
	DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)=DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)+DFstiff(4:6,4:6);
end
for i=2
	DFstiff=Lcrackstiff;
	t=[cos(angle(i)) sin(angle(i)) 0; -sin(angle(i)) cos(angle(i)) 0; 0 0 1];
	T=[t zeros(3,3); zeros(3,3) t];
	DFstiff=T'*DFstiff*T;
	n1=conn(i,1);
	n2=conn(i,2);
	DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)=DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)+DFstiff(1:3,1:3);
	DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)=DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)+DFstiff(1:3,4:6);
	DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)=DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)+DFstiff(4:6,1:3);
	DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)=DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)+DFstiff(4:6,4:6);
end
for i=3:5
	Alpha2(i)=Rho*A*w^2*L(i)^2/(E*A);
	Beta2(i)=Rho*Ip*Mu^2*w^2/(E*A);
	Gamma(i)=(Alpha2(i)/(1-Beta2(i)))^0.5;
	bb(i)=Rho*A*w^2*L(i)^4/(E*I);
	rr(i)=I/(A*L(i)^2);
	ss(i)=E*I/(k*A*G*L(i)^2);
	a1(i)=E*A/L(i)*Gamma(i)*(1-Beta2(i))*cot(Gamma(i));
	a2(i)=-E*A/L(i)*Gamma(i)*(1-Beta2(i))*csc(Gamma(i));
	Phi(i)=(bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	S(i)=sin(Phi(i));
	C(i)=cos(Phi(i));
	W1(i)=E*I/L(i);
	W2(i)=E*I/L(i)^2;
	W3(i)=E*I/L(i)^3;
	Z(i)=Phi(i)-bb(i)*ss(i)/Phi(i);
	if bb(i)*rr(i)*ss(i)<1
		j=1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sinh(Lambda(i));
		C_(i)=cosh(Lambda(i));	
	end
	if bb(i)*rr(i)*ss(i)>1
		j=-1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sin(Lambda(i));
		C_(i)=cos(Lambda(i));	
	end
	Eta(i)=Z(i)/(j*Lambda(i)+bb(i)*ss(i)/Lambda(i));
	Tau(i)=(Lambda(i)+Eta(i)*Phi(i))/(2*Eta(i)*(1-C(i)*C_(i))+(1-j*Eta(i)^2)*S(i)*S_(i));
	d1(i)=W3(i)*bb(i)*Tau(i)*(C(i)*S_(i)+Eta(i)*S(i)*C_(i))/(Lambda(i)*Phi(i));
	d2(i)=W2(i)*Z(i)*Tau(i)*((Phi(i)+j*Eta(i)*Lambda(i))*S(i)*S_(i)-(Lambda(i)-Eta(i)*Phi(i))*(1-C(i)*C_(i)))/(Lambda(i)+Eta(i)*Phi(i));
	d3(i)=W1(i)*Tau(i)*(S(i)*C_(i)-j*Eta(i)*C(i)*S_(i));
	d4(i)=-W3(i)*bb(i)*Tau(i)*(S_(i)+Eta(i)*S(i))/(Lambda(i)*Phi(i));
	d5(i)=W2(i)*Z(i)*Tau(i)*(C_(i)-C(i));
	d6(i)=W1(i)*Tau(i)*(j*Eta(i)*S_(i)-S(i));
	Jm=Jm+int16(fix(Gamma(i)/pi));
	DFstiff=[a1(i) 0 0 a2(i) 0 0; 0 d1(i) d2(i) 0 d4(i) d5(i); 0 d2(i) d3(i) 0 -d5(i) d6(i); a2(i) 0 0 a1(i) 0 0; 0 d4(i) -d5(i) 0 d1(i) -d2(i); 0 d5(i) d6(i) 0 -d2(i) d3(i)];
	t=[cos(angle(i)) sin(angle(i)) 0; -sin(angle(i)) cos(angle(i)) 0; 0 0 1];
	T=[t zeros(3,3); zeros(3,3) t];
	DFstiff=T'*DFstiff*T;
	n1=conn(i,1);
	n2=conn(i,2);
	DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)=DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)+DFstiff(1:3,1:3);
	DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)=DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)+DFstiff(1:3,4:6);
	DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)=DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)+DFstiff(4:6,1:3);
	DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)=DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)+DFstiff(4:6,4:6);
end
for i=6
	DFstiff=Mcrackstiff;
	t=[cos(angle(i)) sin(angle(i)) 0; -sin(angle(i)) cos(angle(i)) 0; 0 0 1];
	T=[t zeros(3,3); zeros(3,3) t];
	DFstiff=T'*DFstiff*T;
	n1=conn(i,1);
	n2=conn(i,2);
	DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)=DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)+DFstiff(1:3,1:3);
	DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)=DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)+DFstiff(1:3,4:6);
	DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)=DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)+DFstiff(4:6,1:3);
	DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)=DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)+DFstiff(4:6,4:6);
end
for i=7:8
	Alpha2(i)=Rho*A*w^2*L(i)^2/(E*A);
	Beta2(i)=Rho*Ip*Mu^2*w^2/(E*A);
	Gamma(i)=(Alpha2(i)/(1-Beta2(i)))^0.5;
	bb(i)=Rho*A*w^2*L(i)^4/(E*I);
	rr(i)=I/(A*L(i)^2);
	ss(i)=E*I/(k*A*G*L(i)^2);
	a1(i)=E*A/L(i)*Gamma(i)*(1-Beta2(i))*cot(Gamma(i));
	a2(i)=-E*A/L(i)*Gamma(i)*(1-Beta2(i))*csc(Gamma(i));
	Phi(i)=(bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	S(i)=sin(Phi(i));
	C(i)=cos(Phi(i));
	W1(i)=E*I/L(i);
	W2(i)=E*I/L(i)^2;
	W3(i)=E*I/L(i)^3;
	Z(i)=Phi(i)-bb(i)*ss(i)/Phi(i);
	if bb(i)*rr(i)*ss(i)<1
		j=1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sinh(Lambda(i));
		C_(i)=cosh(Lambda(i));	
	end
	if bb(i)*rr(i)*ss(i)>1
		j=-1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sin(Lambda(i));
		C_(i)=cos(Lambda(i));	
	end
	Eta(i)=Z(i)/(j*Lambda(i)+bb(i)*ss(i)/Lambda(i));
	Tau(i)=(Lambda(i)+Eta(i)*Phi(i))/(2*Eta(i)*(1-C(i)*C_(i))+(1-j*Eta(i)^2)*S(i)*S_(i));
	d1(i)=W3(i)*bb(i)*Tau(i)*(C(i)*S_(i)+Eta(i)*S(i)*C_(i))/(Lambda(i)*Phi(i));
	d2(i)=W2(i)*Z(i)*Tau(i)*((Phi(i)+j*Eta(i)*Lambda(i))*S(i)*S_(i)-(Lambda(i)-Eta(i)*Phi(i))*(1-C(i)*C_(i)))/(Lambda(i)+Eta(i)*Phi(i));
	d3(i)=W1(i)*Tau(i)*(S(i)*C_(i)-j*Eta(i)*C(i)*S_(i));
	d4(i)=-W3(i)*bb(i)*Tau(i)*(S_(i)+Eta(i)*S(i))/(Lambda(i)*Phi(i));
	d5(i)=W2(i)*Z(i)*Tau(i)*(C_(i)-C(i));
	d6(i)=W1(i)*Tau(i)*(j*Eta(i)*S_(i)-S(i));
	Jm=Jm+int16(fix(Gamma(i)/pi));
	DFstiff=[a1(i) 0 0 a2(i) 0 0; 0 d1(i) d2(i) 0 d4(i) d5(i); 0 d2(i) d3(i) 0 -d5(i) d6(i); a2(i) 0 0 a1(i) 0 0; 0 d4(i) -d5(i) 0 d1(i) -d2(i); 0 d5(i) d6(i) 0 -d2(i) d3(i)];
	t=[cos(angle(i)) sin(angle(i)) 0; -sin(angle(i)) cos(angle(i)) 0; 0 0 1];
	T=[t zeros(3,3); zeros(3,3) t];
	DFstiff=T'*DFstiff*T;
	n1=conn(i,1);
	n2=conn(i,2);
	DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)=DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)+DFstiff(1:3,1:3);
	DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)=DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)+DFstiff(1:3,4:6);
	DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)=DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)+DFstiff(4:6,1:3);
	DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)=DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)+DFstiff(4:6,4:6);
end
DFGstiff(14,14)=DFGstiff(14,14)-M*w^2;
DFGstiff(15,15)=DFGstiff(15,15)-J*w^2;
%Assembling the global dynamic stiffness matrix of the whole system.
for i=1
	if bb(i)*rr(i)*ss(i)<1
		Jm=Jm+int16(fix(Phi(i)/pi))-0.5*(2-sign(d3(i))-sign(d3(i)-d6(i)^2/d3(i)));
	elseif bb(i)*rr(i)*ss(i)==1
		Jm=Jm+int16(fix(Phi(i)/pi))-0.5*(2-sign(d3(i))-sign(d3(i)-d6(i)^2/d3(i)))+1;
	else
		j=-1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		Jm=Jm+int16(fix(Phi(i)/pi))+int16(fix(Lambda(i)/pi))-0.5*(2-sign(d3(i))-sign(d3(i)-d6(i)^2/d3(i)))+1;
	end
end
for i=3:5
	if bb(i)*rr(i)*ss(i)<1
		Jm=Jm+int16(fix(Phi(i)/pi))-0.5*(2-sign(d3(i))-sign(d3(i)-d6(i)^2/d3(i)));
	elseif bb(i)*rr(i)*ss(i)==1
		Jm=Jm+int16(fix(Phi(i)/pi))-0.5*(2-sign(d3(i))-sign(d3(i)-d6(i)^2/d3(i)))+1;
	else
		j=-1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		Jm=Jm+int16(fix(Phi(i)/pi))+int16(fix(Lambda(i)/pi))-0.5*(2-sign(d3(i))-sign(d3(i)-d6(i)^2/d3(i)))+1;
	end
end
for i=7:8
	if bb(i)*rr(i)*ss(i)<1
		Jm=Jm+int16(fix(Phi(i)/pi))-0.5*(2-sign(d3(i))-sign(d3(i)-d6(i)^2/d3(i)));
	elseif bb(i)*rr(i)*ss(i)==1
		Jm=Jm+int16(fix(Phi(i)/pi))-0.5*(2-sign(d3(i))-sign(d3(i)-d6(i)^2/d3(i)))+1;
	else
		j=-1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		Jm=Jm+int16(fix(Phi(i)/pi))+int16(fix(Lambda(i)/pi))-0.5*(2-sign(d3(i))-sign(d3(i)-d6(i)^2/d3(i)))+1;
	end
end
if isempty(sup1)&&isempty(sup2) %In the case of no supports or restraints.
    DFGstiff=genre(DFGstiff); %Perfom Gaussian elimination using the function file 'genre.m'.
    diagonal=diag(DFGstiff);
    for i=1:length(diagonal);
        if diagonal(i)<0
	        ng=ng+1; %Sign counts of the diagonal terms.
	    end
    end
    WW=ng+Jm; %Summation of terms of the Wittrick-Williams algorithm.
else
    DFGstiff([sup2],:)=[];
    DFGstiff(:,[sup2])=[];
	DFGstiff([sup1],:)=[];
	DFGstiff(:,[sup1])=[];%Delete the rows and columns in the stiffness matrix, corresponding to the suppressed degrees of freedom.
    DFGstiff=genre(DFGstiff);
    diagonal=diag(DFGstiff);
    for i=1:length(diagonal);
        if diagonal(i)<0
	        ng=ng+1;
	    end
    end
    WW=ng+Jm;
end
end