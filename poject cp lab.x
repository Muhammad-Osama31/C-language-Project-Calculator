#include<stdio.h>
#include<conio.h>
void main()
{
	printf("\n-----------------------------------------------------------");
	printf("\nThis is physics Calulator which formula your want to solve ");
	printf("\n-----------------------------------------------------------");
	printf("\n Enter which formula you want\n");
	printf("\n-------------------------------");
	printf("\n For  Velocity plz Enter :-> 1  \n");
	printf("\n For Acceleration plz Enter :->2 \n");
	printf("\n For Weight Formula plz Enter:->3 \n");
	printf("\n For Newtons Second Law plz Enter :->4 \n");
	printf("\n For Centripital Accerleration plz Enter :->5 \n");
	printf("\n For Momentum plz Enter value  :->6 \n");
	printf("\n For Efficency plz Enter value :-> 7 \n");
	printf("\n For power plz Enter value :->8 \n");
	printf("\n For Kinetic Energy plz Enter :->9 \n");
	printf("\n For Average  Force formua plz Enter:->10 \n");
	printf("\n For Amplitude plz Enter :->11 \n");
	printf("\n For Acceleration due to gravity plz Enter :->12 \n");
	printf("\n For Archimedes Principle Formula plz Enter :->13 \n");
	printf("\n For Enropy is plz Enter value :->14 \n");
	printf("\n For  Efficency of Carnot Engine plz Enter value :->15 \n");
	printf("\n For Coulomb's law is plz Enter value :->16 \n");
	printf("\n For Electron volt plz Enter value:->17 \n");
	printf("\n For Ohm's Law plz Enter value :->18 \n");
	printf("\n For  Series Combination plz Enter value :->19 \n");
	printf("\n For Resistance of Series plz Enter value :->20 \n");
	printf("\n For  Diameter of circle plz Enter value :->21 \n");
	printf("\n For Circumference of circle plz Enter value :->22 \n");
	printf("\n For Area of circle plz Enter value :->23 \n");
	printf("\n For Curved surface area of a cone  plz Enter value:->24 \n");
	printf("\n For Current Density Formula plz Enter value:->25 \n");
	printf("\n For Area of Sphere Formula plz Enter value :-26 \n");
	printf("\n For Speed Distance Time Formula Plz Enter Value :->27 \n");
	printf("\n For  DC voltage drop formula Plz Enter Value :->28 \n");
	printf("\n For The Gravity Formula  Plz Enter Value :->29 \n");
	printf("\n For Latent Heat Formula  Plz Enter Value :->30 \n");
	printf("\n For The Youngs modulus formula  Plz Enter Value :->31 \n");
	printf("\n For The UpWard Tension formula  Plz Enter Value :->32 \n");
	printf("\n For The Downward Tension Formula  Plz Enter Value :->33 \n");
	printf("\n For wave formula Plz Enter Value :->34 \n");
	printf("\n For The wavelenght formula Plz Enter Value :->35 \n");
	printf("\n For The Moment of inertia Plz Enter Value :->36 \n");
	printf("\n For uniform circular motionz plz Enter Value :->37 \n"); 
	printf("\n For The Frequency plz Enter Value :->38 \n");
	printf("\n For The Timeperiod plz Enter Value :->39 \n");
	printf("\n For The Mass Energy plz Enter Value :->40 \n");
	printf("\n For The Density Formula plz Enter Value :->41 \n");
	printf("\n For The Electric Charge Formula plz Enter Value :->42 \n");
	printf("\n For The Torque  Formula plz Enter Value :->43 \n");
	printf("\n For Gravitaional potential Energy  plz Enter Value :->44 \n");
	printf("\n For The Work done plz Enter Value :->45 \n");
	printf("\n For The Specific Heat formula is plz Enter Value :->46 \n");
	printf("\n For The Capacitance of capacitor  plz Enter Value :->47 \n");
	printf("\n For Convertions of Celsius into Fahrenheit plz Enter Value :->48 \n");
	printf("\n For Convertions of Fahrenheit into Celsius  plz Enter Value :->49 \n");
	printf("\n Convertions of Celsius into Kelvin plz Enter Value :->50  \n\n\n");
	int i;
	scanf_s("%d", &i);
	
		if (i == 1)
		{
			float s, t, v;
			printf("\nFor velocity please Enter value of s and t v=s/t: \n");
			scanf_s("%f%f", &s, &t);
			v = s / t;
			printf("\nThe velocity is %f", v);
			_getch();
			main();
		}
		if (i == 2){
			float a, v, t;
			printf("\nFor accerlation please Enter value of v and t a=v/t \n");
			scanf_s("%f%f", &v, &t);
			a = v / t;
			printf("\nThe Acceration is %f ", a);
			_getch();
			main();
		}
		if (i == 3){

			float w, m, g;
			printf("\nFor wight please Enter value of m and g w=m*g \n");
			scanf_s("%f%f", &m, &g);
			w = m*g;
			printf("\nThe Weight formula is %f", w);
			_getch();
			main();
		}
		if (i == 4)
		{
			float f, m, a;
			printf("\nFor Newtons Second Law please Enter value of m and a  f=m*a \n");
			scanf_s("%f%f", &m, &a);
			f = m*a;
			printf("\n The Newtons Seond law is %f", f);
			_getch();
			main();
		}

		if (i == 5)
		{
			float  a, v, r;
			printf("\nFor Centripital Accerleration here v=velocity andr=radius a=v^2/r \t\tplz Enter value of v and r \n");
			scanf_s("%f%f", &v, &r);
			a = v*v / r;
			printf("\n The Centipital Accerleration is %f", a);
			_getch();
			main();
		}
		if (i == 6){
			float p, m, v;
			printf("\nFor Momentum  here the m=mmomentum and v=velocity Formula is p=m*v \t\tplz Enter value of v and m\n");
			scanf_s("%f%f", &m, &v);
			p = m*v;
			printf("\n The Momentum is %f", p);
			_getch();
			main();
		}
		if (i == 7){
			float Q, e, W;
			printf("\nFor Efficency here of W=workdont and Q=Heat Formula is e=W/Q \t\t\tplz Enter value of W and Q \n");
			scanf_s("%f%f", &W, &Q);
			e = W / Q;
			printf("\n The Efficency is %f", e);
			_getch();
			main();
		}
		if (i == 8)
		{
			float p, w, t;
			printf("\nFor power here W=work and t=time Formula is p=w/t plz Enter value of W and T\n");
			scanf_s("%f%f", &w, &t);
			p = w / t;
			printf("\n The Power is %f", p);
			_getch();
			main();
		}
		if (i == 9)
		{
			float k, m, v;
			printf("\nFor Kinetic Energy of m=mass and v=velocity \n k=1/2*m*v^2 \nplz Enter value m and v  \n");
			scanf_s("%f%f", &m, &v);
			k = 0.5*m*v*v;
			printf("\n The kinetic Energy is %f", k);
			_getch();
			main();
		}
		if (i == 10)
		{
			float F, m, vi, vf, dt;
			printf("\n For Average  Force formua \n m=mass ,vi=initital momentum  ,vf=finial momentum and dt=change in time \t F=m(vf-vi)/dt plz Enter value of m,vf,vi,dt \n");
			scanf_s("%f%f %f %f", &m, &vf, &vi, &dt);
			F = m*(vf - vi) / dt;
			printf("\nThe Average Force Formula is %f ", F);
			_getch();
			main();
		}
		if (i == 11)
		{
			float A, D, F;
			printf("\nFor Amplitude plz Enter value of D=distance and F=frequency  A=D/F \n");
			scanf_s("%f%f", &D, &F);
			A = D / F;
			printf("\nThe Amplitude is %f", A);
			_getch();
			main();
		}
		if (i == 12)
		{
			float g, G, M, r, h;
			printf("\n For Acceleration due to gravity  G= is the universal gravitational,\n r is the radius,h is the height, m= is the mass , g= GM/(r+h)^2 Plz Enter value of G,M,r,h \n\n");
			scanf_s("\n %f%f%f%f", &G, &M, &r, &h);
			g = G*M / (r + h) * 2;
			printf("\nThe Acceleration due to gravity is %f", g);
			_getch();
			main();
		}
		if (i == 13)
		{
			float F, p, g, V;
			printf("\n For Archimedes Principle Formula Buoyant force of a given body = F,\n Volume of the displaced fluid = v ,acceleration due to gravity = g F=pgv \nplz Enter value of p,v,g \t\n");
			scanf_s("\n %f%f%f", &p, &g, &V);
			F = p*g*V;
			printf("\n Archimedes Principle is %f ", F);
			_getch();
			main();
		}
		if (i == 14)
		{
			float dQ, T, dS;
			printf("\n For Enropy is where dt=temperature ,dq=heat transfered in system,  dS=dQ/T plz Enter value dQ and T\n");
			scanf_s("%f%f", &dQ, &T);
			dS = dQ / T;
			printf("\nThe Entropy is %f", dS);
			_getch();
			main();
		}
		if (i == 15)
		{
			float e, T1, T2;
			printf("\nFor Efficency of Carnot Engine is t2=low temperature and t1=high temperature plz Enter value  e=1-T2/T1*100\n plz Enter value of T2 and T1\n");
			scanf_s("%f%f", &T2, &T1);
			e = 1 - T2 / T1 * 100;
			printf("\nThe Efficiency of Carnot Engine is %f", e);
			_getch();
			main();
		}
		if (i == 16)
		{
			float q1, q2, r, k, F;
			printf("\nFor Coulomb's law\nq1=is heat observed by system ,\nq2=is heat rejected by the system ,\nk=constant of proportionality F=k q1*q2/r^2 \nplz Enter value of k,q1,q2,r\n");
			scanf_s("%f%f%f%f", &k, &q1, &q2, &r);
			F = k*(q1*q2) / r*r;
			printf("\n The Coulomb's Law is %f\n", F);
			_getch();
			main();
		}
		if (i == 17)
		{
			float m, v, E;
			printf("\n For Electron volt plz Enter value of M=mass of electron V=potentital difference 1/2*m*v plz Enter value of M and V\n");
			scanf_s("%f%f", &m, &v);
			E = 0.5*m*v;
			printf("\n The Electron volt is %f", E);
			_getch();
			main();
		}
		if (i == 18)
		{
			float I, V, R;
			printf("\n For Ohm's Law plz Enter value where V=voltage and R=resistance I=V/R  \t\tplz Enter value of V and R\n");
			scanf_s("%f%f", &V, &R);
			I = V / R;
			printf("\nThe Ohm's Law is %f ", I);
			_getch();
			main();
		}
		if (i == 19)
		{
			float c, c1, c2, c3;
			printf("\n The Series Combination is 1/c = 1/c1+1/c2+1/c3\n where C=is the Capacitance of Each Plate\n plz Enter value of c1, c2, c3\n");
			scanf_s("%f%f%f", &c1, &c2, &c3);
			c = 1 / c1 + 1 / c2 + 1 / c3;
			printf("\n The Series Combination is 1/c is %f ", c);
			_getch();
			main();
		}
		if (i == 20)
		{
			float R, R1, R2, R3;
			printf("\n The Resistance is 1/R = 1/R1+1/R2+1/R3 \nwhere R=Resistance of Each plate \nplz Enter value of R1, R2, R3 \n");
			scanf_s("%f%f%f", &R1, &R2, &R3);
			R = 1 / R1 + 1 / R2 + 1 / R3;
			printf("\n The Resistance is 1/R is %f", R);
			_getch();
			main();
		}
		if (i == 21)
		{
			float r, d;
			printf("\n The Diameter of circle is D=2*r \nplz Enter value of r \n");
			scanf_s("%f", &r);
			d = r * 2;
			printf("\n The Diameter of circle is D=%f", d);
			_getch();
			main();
		}
		if (i == 22)
		{
			float pi, r, C;
			printf("\nCircumference of circle = 2 *pi*r\nplz Enter value of pi and r\n");
			scanf_s("%f%f", &pi,&r);
			C = 2 * pi*r;
			printf("\nCircumference of circle is 2 *pi*r= is %f ", C);
			_getch();
			main();
		}
		if (i == 23)
		{
			float A, r, pi;
			printf("\nArea of circle = 2*pi*r2 \nplz Enter value of pi and r \n");
			scanf_s("%f%f",&pi,&r);
			A = 2 * pi*r*r;
			printf("\nArea of circle  2*pi*r2 is=%f ", A);
			_getch();
			main();
		}
		if (i == 24)
		{
			float A, r, l, pi;
			printf("\n Curved surface area of a cone is pi*r*l \nEnter value of pi , r and l\n");
			scanf_s("%f%f%f",&pi, &r, &l);
			A = pi*r*l;
			printf("\nCurved surface area of a cone is %f", A);
			_getch();
			main();
		}
		if (i == 25)
		{
			float J, I, A;
			printf("\nCurrent Density Formula is expressed as,J = I / A \nEnter value of I and A\n");
			scanf_s("%f%f", &I, &A);
			J = I / A;
			printf("\nCurrent Density is %f", J);
			_getch();
			main();
		}
		if (i == 26)
		{
			float pi, r, A;
			printf("\nThe Area of Sphere is expressed as,3/4*pi*r^2\nEnter value of pi and r\n");
			scanf_s("%f%f", &pi, &r);
			A = 0.75* pi*r*r;
			printf("\nThe Area of Sphere is %f", A);
			_getch();
			main();
		}
		if (i == 27)
		{
			float s, x, t;
			printf("\nSpeed Distance Time Formula is x = Speed in m/s,\nd = Distance traveled in m,\nt = time ,taken in speed.distance = X*t \nEnter value of x and t\n");
			scanf_s("%f%f", &x, &t);
			s = x*t;
			printf("\nSpeed Distance Time is %f", s);
			_getch();
			main();
		}
		if (i == 28)
		{
			float V, L, I, T;
			printf("\nDC voltage drop formula isV=L*I/T\nWhere I = current through the circuit,\nL = length of the circuit \nEnter value of L,I,T\n");
			scanf_s("%f%f%f", &L, &I, &T);
			V = L*I / T;
			printf("\nThe DC voltage drop formula is %f ", V);
			_getch();
			main();
		}
		if (i == 29)
		{
			float  m1, m2, r, F, G,A,E;
			printf("\nThe Gravity Formula is F=G*m1*m2/r^2 G is a constant equal to 6.67 * 10-11 N-m2/kg2,\nEnter value of G,\nm1 is the mass of the body 1,\nm2 is the mass of body 2,\nr is the radius\ndistance amid the two bodies\n Please Enter value ofG,m1,m2,r\n");
			scanf_s("%f%f%f%f", &G,&m1,&m2,&r);
			A = G*m1*m2 ;
			F = r*r;
			E = A / F;
			printf("\nThe Gravity Formula is %f \n", E);
			_getch();
			main();
		}
		if (i == 30)
		{
			float Q, m, T, c;
			printf("\n Latent Heat Formula is c=Q/m*T Where,\nm = mass of the body,\nC = specific heat, \nT= is the temperature difference \nplz Enter value of Q,m,T\n");
			scanf_s("%f%f%f", &Q, &m, &T);
			c = Q / m*T;
			printf("\nHeat Formula is %f\n", c);
			_getch();
			main();
		}
		if (i == 31)
		{
			float E, e, fi;
			printf("\nThe Youngs modulus formula is E=fi/e\nWhereE = Young’s modulus,\ns = stress,e strain \nplz Enter value of fi and e\n");
			scanf_s("%f%f", &fi, &e);
			E = fi / e;
			printf("\nThe Youngs modulus is %f\n", E);
			_getch();
			main();
		}
		if (i == 32)
		{
			float T, W, M, a;
			printf("\n The UpWard Tension Formula is T = W+M*a \nWhere The Weight of the body = W,\nmass of the body = m,\nacceleration of the moving body = a \nplz Enter W,M,a\n");
			scanf_s("%f%f%f", &W, &M, &a);
			T =W + M*a;
			printf("\nThe UpWard Tension is %f\n", T);
			_getch();
			main();
		}
		if (i == 33)
		{
			float T, W, M, a;
			printf("\n The Downward Tension Formula is T = W-M*a \n3Where,The Weight of the body = W,\nmass of the body = m,\nacceleration of the moving body = a \nplz Enter W,M,a\n");
			scanf_s("%f%f%f", &W, &M, &a);
			T = W - M*a;
			printf("\nThe DownWard Tension is %f\n", T);
			_getch();
			main();
		}

		if (i == 34)
		{
			float V, f, l;
			printf("The wave formula is V = f *lamda Where,\nv = velocity of the wave,\nf = frequency of the wave,\nlamda = wavelength\nPlz Enter value of f,lamda\n");
			scanf_s("%f%f", &f, &l);
			V = f* l;
			printf("\n The Wave formula is %f", V);
			_getch();
			main();
		}
		if (i == 35)
		{
			float V, f, l;
			printf("The wavelenght formula is V = f /lamda Where,\nv = velocity of the wave,\nf = frequency of the wave,\nlamda = wavelength\nPlz Enter value of f,lamda\n");
			scanf_s("%f%f", &f, &l);
			V = f / l;
			printf("\n The Wave formula is %f", V);
			_getch();
			main();
		}
		if (i == 36)
		{
			float I, M, R;
			printf("\nMoment of inertia formula is given by I = MR^2 Where,\nR is the Distance between the axis and rotation in m,\nM is the Mass of the object in Kg,\nI = Moment of Inertia in Kgm^2,\nEnter value of M and R\n");
			scanf_s("%f%f", &M, &R);
			I = M*R*R;
			printf("\nMoment of inertia  is %f", I);
			_getch();
			main();
		}
		if (i == 37)
		{
			float V, T, pi , r;
			printf("\nuniform circular motion is V = 2*pi*r/ T \nplz Enter value of pi,r,T\n");
			scanf_s("%f%f%f",&pi, &r, &T);
			V = 2 * pi*r / T;
			printf("\nuniform circular motion is %f", V);
			_getch();
			main();
		}
		if (i == 38)
		{
			float F, T;
			printf("\nThe Frequency is f = 1/T \nplz Enter value of T\n");
			scanf_s("%f", &T);
			F = 1 / T;
			printf("\nThe Frequency is %f\n", F);
			_getch();
			main();
		}
		if (i == 39)
		{
			float F, T;
			printf("\nThe Timeperiod is T = 1 / F \nplz Enter value of F\n");
			scanf_s("%f", &F);
			T = 1 / F;
			printf("\nThe Timeperiod is %f\n", T);
			_getch();
			main();
		}
		if (i == 40)
		{
			float E, M,C;
			printf("\nThe Mass Energy is E=M*C^2 \nplz Enter value of M and C\n");
			scanf_s("%f%f", &M,&C);
			E=M*C*C;
			printf("\nThe Mass Energy is %f\n", E);
			_getch();
			main();
		}
		if (i == 41)
		{
			float P,m,v;
			printf("\nThe Density Formula is P = m / v \nplz Enter value of m and v\n");
			scanf_s("%f%f", &m, &v);
			P = m / v;
			printf("\nThe Density is %f\n", P);
			_getch();
			main();
		}
		if (i == 42)
		{
			float Q,I,T;
			printf("\nThe Electric Charge Formula is Q = I*T \nplz Enter value of I and T\n");
			scanf_s("%f%f", &I,&T);
			Q = I*T;
			printf("\nThe Electric Charge is %f\n", Q);
			_getch();
			main();
		}
		if (i == 43)
		{
			float F,r, T;
			printf("\nThe Torque  Formula is T = r*F; \nplz Enter value of r and F\n");
			scanf_s("%f%f", &r, &F);
			T = r*F;
			printf("\nThe Torque is %f\n", T);
			_getch();
			main();
		}if (i == 44)
		{
			float E, m, g, h;
			printf("\nGravitaional potential Energy is E = m*g*h \nplz Enter value of m,g,h \n");
			scanf_s("%f%f%f", &m, &g,&h);
			E = m*g*h;
			printf("\nGravitaional potential Energy is %f\n", E);
			_getch();
			main();
		}
		if (i == 45)
		{
			float W, f, d;
			printf("\nThe Work done Formula is W= f*d \nplz Enter value of f,d \n");
			scanf_s("%f%f", &f, &d);
			W= f*d;
			printf("\nThe Work done is %f\n", W);
			_getch();
			main();
		}
		if (i == 46)
		{
			float c, dQ, m,dT,q;
			printf("\nThe Specific Heat formula is c=dQ/m*dT\nplz Enter value of dQ,m,dT \n");
			scanf_s("%f%f%f", &dQ, &m,&dT);
			q = m*dT;
			c = dQ / q;
			printf("\nThe Specific Heat is %f\n", c);
			_getch();
			main();
		}
		if (i == 47)
		{
			float C,q,v;
			printf("\nThe Capacitance of capacitor Formula is C=q/v; \nplz Enter value of q,v \n");
			scanf_s("%f%f", &q, &v);
			C = q / v;
			printf("\nCapacitance of capacitoris %f\n", C);
			_getch();
			main();
		}if (i == 48)
		{
			float  T, F;
			printf("\nConvertions of Celsius into Fahrenheit Formula is F=2/5*T+32 \nplz Enter value of T \n");
			scanf_s("%f", &T);
			F = 0.4*T+32;
			printf("\nConvertions of Celsius into Fahrenheit is %f\n", F);
			_getch();
			main();
		}
		if (i == 49)
		{
			float  T, C;
			printf("\nConvertions of Fahrenheit into Celsius Formula is C=5/9*T-32 \nplz Enter value of T \n");
			scanf_s("%f", &T);
			C = 0.5555 * T - 32;
			printf("\nConvertions of Fahrenheit into Celsius %f\n", C);
			_getch();
			main();
		}
		if (i == 50)
		{
			float  T, K;
			printf("\nConvertions of Celsius into Kelvin Formula is C=T+273 \nplz Enter value of T \n");
			scanf_s("%f", &T);
			K =  T +273;
			printf("\nConvertions of Celsius into Kelvin is %f\n", K);
			_getch();
			main();
		}

	_getch();

}