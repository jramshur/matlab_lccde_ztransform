%% LCCDE and Z-Transforms with MATLAB
% Purporse Statment: The purpose of this document is to
% demonstrate how use MATLAB in working with LCCDE's,
% z-transforms, and system functions. MATLAB will be used to
% 
% # solve LCCDE's
% # determine zeros, poles, and zero-pole diagrams
% # determine impulse and step responses.
%
% Please consult any function documentation for details of usage. The methods 
% discussed here are not an exhastive list of methods.

%% Linear Constanst-Coeficient Difference Equations (LCCDE)
%
% An important class of LTI systems consist of those systems for which the input x[n]
% and the output y[n] satisfy and Nth-order linear constant-coefficient difference equation
% of the form
% 
% $$ y[n] + a_1y[n-1] + ... + a_Ny[n-N] = b_0x[n] + b_1x[n-1] + ... + b_Mx[n-M] $$
%
% or
%
% $$ \sum_{k=0}^{N}a_ky[n-k] = \sum_{k=0}^{M}b_kx[n-k]. $$
% 
% This can also be arranged in the recursive form
%
% $$ y[n]=-\sum_{k=1}^{N}(\frac{a_k}{a_0})y[n-k] + \sum_{k=0}^{M}(\frac{b_k}{a_0})x[n-k]. $$
% 
%
% Solutions to LLCDE's can be obtained using various methods. We will discuss
% methods involvings recursion and z-transforms. Solving LCCDE's using
% the recursion method requires prior knoledge of inintal conditions
% and of the function input sequence x[n].
%
% *Example 1: Recursion*
% 
% (Oppenheim, Problem 2.5)
% 
% Solve the following LLCDE by recursion.
%
% $$ y[n]-5y[n-1]+6y[n-2]=2x[n-1] $$
%
% Given that x[n]=u[n], y[n]=0 for n<0, find y[n]. The solution to this will be
% the system step response. The book's solution for a step response is 
% $$ h[n]=u[n]-4(2)^nu[n]+3(3)^nu[n] $$ .

% solve recursively
a0=1; % a_0 term
a=[-5 6]./a0; % 'a_k' terms, k=1:N
b=[0 2 0]./a0;  % 'b_k' terms, k=0:M
p=length(a);    % number of a_k terms, order of system
q=length(b);    % number of b_m terms
N=20;
x=[zeros(1,p) ones(1,N+1)]; % input sequence (step sequence)
y=zeros(1,N+p+1); % preallocate memory for output sequence
for n=p+1:length(y)    
    y(n)= -sum(a.*y(n-1:-1:n-p)) + sum(b.*x(n:-1:n-q+1)); %compute y[n]
end
y=y(p+1:end); % keep values from n = 0 to N
figure;
stem(subplot(211),[0:20],y) % plot
title('Step Response (recursive solution)'); 
xlim([0 21]); xlabel('n (samples)'); ylabel('Amplitude');

% Compare with book solution
n=0:20;
h = 1 - 4*2.^n + 3*3.^n;    % compute solution
stem(subplot(212),[0:20],h) % plot
title('Step Response (book solution)'); 
xlim([0 21]); xlabel('n (samples)'); ylabel('Amplitude');
%%
% Given that x[n]=delta[n], y[n]=0 for n<0, find y[n]. The solution will be
% the system impulse response. The book's solution for an impulse response is
% $$ h[n]=-2(2)^nu[n]+2(3)^nu[n] $$ .

% solve recursively
x=[zeros(1,p) 1 zeros(1,N)]; % input sequence (impulse sequence)
y=zeros(1,N+p+1); % preallocate memory for output sequence
for n=p+1:length(y)    
    y(n)= -sum(a.*y(n-1:-1:n-p)) + sum(b.*x(n:-1:n-q+1)); %compute y[n]
end
y=y(p+1:end); % keep values from n = 0 to N
figure;
stem(subplot(211),[0:20],y) % plot
title('Impulse Response (recursive solution)'); 
xlim([0 21]); xlabel('n (samples)'); ylabel('Amplitude');

% Compare with book solution
n=0:20;
h = -2*2.^n + 2*3.^n;    % compute solution
stem(subplot(212),[0:20],h) % plot
title('Impulse Response (book solution)'); 
xlim([0 21]); xlabel('n (samples)'); ylabel('Amplitude');


%% z-Transform
%
% The z-transform allows a more direct and/or closed form solution to a 
% difference equation compared to recusion. Also, the z-transform allows us
% to extract the zeros and poles of a system.
%
% The z-transform of a sequence x[n] is defined as
%
% $$ X(z)=\sum_{n=-\infty}^{\infty}x[n]z^{-n}=Z[x[n]]. $$
%
% Recall the LCCDE form
%
% $$ y[n]=-\sum_{k=1}^{N}(\frac{a_k}{a_0})y[n-k] + \sum_{k=0}^{M}(\frac{b_k}{a_0})x[n-k]. $$
%
% Taking the z-tranform of both sides yields
%
% $$ Y(z)=-\sum_{k=1}^{N}(\frac{a_k}{a_0})z^{-k}Y(z) + \sum_{k=0}^{M}(\frac{b_k}{a_0})z^{-k}X(z). $$
%
% Solving Y(z) in terms of X(z) yields
%
% $$ Y(z)=\frac{\sum_{k=0}^{M}{b_k}z^{-k}}{\sum_{k=0}^{N}{a_k}z^{-k}}X(z)=H(z)X(z) $$
%
% where H(z) represents the system function of the LTI system.
%
% $$ H(z)=\frac{Y(z)}{X(z)} $$ 
%
%% z-Transforms in MATLAB
%
% MATLAB includes several functions that can aid in performing z-transforms (ZT) 
% of a function or sequence and/or it's inverse (iZT). Three methods and their respective
% functions are discussed below.
%
% <html><h3>Symbolic Toolbox Method</h3></html>
% 
% The _ztrans()_ and _iztrans()_ functions perform the
% z-transform and the inverse z-transform.
%
% *Example 2: Sympolic ZT*
% 
% Find ZT of x[n]=a^n u[n]
%
syms a n x X % define symbolic variables
x=a^n;       % define sequence x[n] (Table 3.1)
X=ztrans(x)  % perform ZT

%%
% *Example 3: Sympolic ZT*
%
% Find ZT of x[n]=cos[wn]
%
syms w n x X % define symbolic variables
x=cos(w*n);  % define sequence x[n] (Table 3.1)
X=ztrans(x)  % perform ZT

%%
% *Example 4: Sympolic iZT*
% (Oppenhiem, Ex. 3.10) 
%
% Find iZT of 
%
% $$ X(z)=\frac{1+2z^{-1}+z^{-2}}{1-\frac{3}{2}z^{-1}+\frac{1}{2}z^{-2}} $$
%
syms z x X   % define symbolic variables
X=(z^2+2*z+1)/(z^2-3*z/2+1/2) %define transfer function X(z)
x=iztrans(X) % perform inverse ZT

%%
% <html><h3> 2. Power Series Method </h3></html>
%
% The *deconv* function can be used to perform the long
% division required in power series method.
%
% For the given z-transform
%
% $$ X(z)=\frac{b_0+b_1z^{-1}+...+b_Mz^{-M}}{a_0+a_1z^{-1}+...+a_Mz^{-M}} $$
%
% the matlab command is *[q,r]=deconv(b,a)* where b and a vectors
% containing the coefficients of the numerator and denominator
% respectively.
%
% *Example 5: iZT by deconvolution*
% 
% Given X(z) find x[n] by deconvolution
%
% $$ X(z)=\frac{1+2z^{z-1}+z^{-2}}{1-z^{z-1}+0.4z^{-2}} $$

b=[1 2 1];
a=[1 -1 .4];
n=20;
b=[b zeros(1,n)]; % make b longer
[x,r]=deconv(b,a);
% plot
figure;
stem(subplot(211),[0:length(x)-1],x) % plot
title('Power Series Expansion Solution of x[n]'); xlim([0 21])
ylabel('Amplitude')
% Compare impulse response
subplot(212);
impz(b,a); % plot impulse response
title('Impulse Response x[n]'); xlim([0 21])

%% 
% *3. Partial fraction expansion (PFE) method*
% 
% $$ H(z)=\frac{B(z)}{A(z)}=\frac{b_0+b_1z^{-1}+...+b_Mz^{-M}}{a_0+a_1z^{-1}+...+a_Mz^{-M}} $$
% 
% We can use the MATLAB *residuez* function to convert between
% partial fraction expansion and polynomial coefficients. Partial fraction expansion
% may yield a z-function that is then solvable using the inspection method. 
% If there are no multiple roots, then
%
% $$ \frac{B(z)}{A(z)}=\frac{r_1}{1-p_1z^{-1}} +\frac{r_2}{1-p_2z^{-1}}+...+\frac{r_n}{1-p_nz^{-1}}+k_1+k_2z^{-1}+...+k_{m-n-1}z^{-m-n} $$
%
% *Example 8: PFE*
%
% Find the partial fraction expansion of X(z) that is given by
%
% $$ X(z)=\frac{1}{(1-\frac{1}{4}z^{-1})(1-\frac{1}{2}z^{-1})} $$
% 
b=1; % numerator coefficients
a=conv([1 -1/4],[1 -1/2]); % multiplying polynomials to get denom. coef.
[r,p,k]=residuez(b,a) % compute residues (r), poles (p), and direct terms (k)

%%
% This represents the partial fraction expansion
%
% $$ X(z)=\frac{2}{1-\frac{1}{2}z^{-1}}+\frac{-1}{1-\frac{1}{4}z^{-1}} $$

%%
% *Example 7: PFE*
% (Oppenheim, Ex. 3.10)
%
% Find the partial fraction expansion of X(z) that is given by
%
% $$ X(z)=\frac{1+2z^{-1}+z^{-2}}{1-\frac{3}{2}z^{-1}+\frac{1}{2}z^{-2}} $$
% 
b=[1 2 1]; % numerator coefficients
a=[1 -3/2 1/2]; % denominator coefficients
[r,p,k]=residuez(b,a) % compute residues (r), poles (p), and direct terms (k)

%%
% This represents the fraction expansion
%
% $$ X(z)=\frac{8}{1-z^{-1}}+\frac{-9}{1-\frac{1}{2}z^{-1}}+2 $$
%

%% Pole-Zero Diagram
%
% Matlab contains several functions that may help when working with zeros and
% poles. The *zplane* and *pzmap* functions can be used display pole-zero (PZ)
% diagram of z-functions. The poles and zeros can be determined from 
% rational polynomials using *roots*, *poles*, and *zero* functions.
%
% *Example 8: Poles, Zeros, and PZ map* 
%
% Given H(z) below find the zeros, poles, and plot the PZ map. We
% will use three methods to find the poles and zeros, and use
% three methods to plot the PZ map.
%
% $$ H(z)=\frac{1-2z^{-1}}{1-\frac{2}{3}z^{-1}} $$
%

a=[1 -2/3]; % denominator coefs.
b=[1 -2];   % numerator coefs.

% ZP Method 1 - using 'roots' and 'zplane' functions
p1=roots(a);    % get poles of rational function
z1=roots(b);    % get zeros of rational function
disp(['Method 1: poles=[' num2str(p1) '] , zeros=[' num2str(z1) ']'])

% ZP Method 2 - using 'tf', 'pole', 'zero', and 'pzplot' functions
H=tf(b,a);      % generate trasfer function based on rational polynomial
p2=pole(H);     % get poles of transfer function
z2=zero(H);     % get zeros of transfer function
disp(['Method 2: poles=[' num2str(p2) '] , zeros=[' num2str(z2) ']'])
H               % print transfer funtion to command line 

% ZP Method 3 - using 'tf2zpk'
[z3,p3,k3] = tf2zpk(b,a);   %Convert transfer function filter parameters to zero-pole-gain form
disp(['Method 3: poles=[' num2str(p3) '] , zeros=[' num2str(z3) '], gain=' num2str(k3)])

% Plot Method 1: using zplane(b,a). NOTE: a and b must be row vectors.
figure; zplane(b,a);  
title('PZ Map using zplane(b,a)'); 
xlim([-2 2]); ylim([-2 2]);
xlabel('Real Axis'); ylabel('Imaginary Axis')

% Plot Method 2: using zplane(z,p). NOTE: z and p must be column vectors.
figure; zplane(z1,p1);  
title('PZ Map using zplane(z,p)'); 
xlim([-2 2]); ylim([-2 2]);
xlabel('Real Axis'); ylabel('Imaginary Axis')

% Plot Method 3: using pzmap
figure; pzmap(H);     
title('PZ Map using pzmap');
xlim([-2 2]); ylim([-2 2]);
xlabel('Real Axis'); ylabel('Imaginary Axis')

%% Problem 1: Convolution
% (Ingle, Ex. 2.10)
%
% Consider the finite duration input sequence x[n] given by
%
% $$ x[n]=u[n]-u[n-10] $$
% 
% and the infinite duration impulse response h[n] is given by
%
% $$ h[n]=(0.9)^nu[n]. $$
%
% Determine y[n]=x[n]*h[n]
%

n=-5:50;            % define range of n
x = stepseq(0,-5,50) - stepseq(10,-5,50); % generate x[n]
h=.9.^(0:50);       % generate h[n]
y_conv=conv(x,h);   % convolve x and h
y_conv=y_conv(1:length(x)); % keep samples of interest
y_filt=filter(1,[1 -.9],x); % compare with filter function

% plot convoltion
figure; subplot(211);
stem(n, y_conv); % plot output sequence
xlim([-5 50]); 
title('Output Sequence y[n] (convolution)'); 
xlabel('n (samples)'); ylabel('Amplitude');

% plot filter comparison
subplot(212);
stem(n, y_filt); % plot output sequence
xlim([-5 50]); 
title('Output Sequence y[n] (filter fxn)'); 
xlabel('n (samples)'); ylabel('Amplitude');
%% Problem 2: iZT
%
% Determine the iZT of X(z) so that the resulting sequence is
% causual and contains no complex numbers where X(z) is given as
%
% $$ X(z)=\frac{1+0.4\sqrt{2}z^{-1}}{1-0.8\sqrt{2}z^{-1}+0.64z^{-2}}. $$
%
b=[1 0.4*sqrt(2)];
a=[1 -0.8*sqrt(2) 0.64];
[r,p,c] = residuez(b,a)
Mp=abs(p')      % pole magnitude
Ap=angle(p')/pi % pole angles in pi units

%%
% From the above we get
%
% $$ X(z)=\frac{0.5+j}{1-\left | 0.8 \right | e^{-j\frac{\pi}{4}}z^{-1}} + \frac{0.5-j}{1-\left | 0.8 \right | e^{j\frac{\pi}{4}}z^{-1}} , \left | z \right | > 0.8 $$
%
% and from Table 3.1 we have
% 
% $$ x[n]=(0.5+j) \left | 0.8 \right | ^ne^{-j\frac{\pi}{4}n}u[n]] + (0.5-j) \left | 0.8 \right | ^ne^{j\frac{\pi}{4}n}u[n]  $$
%
% $$ x[n]=\left | 0.8\right |^n  \left [ 0.5 \left \{ e^{-j\frac{\pi}{4}n} + e^{j\frac{\pi}{4}n} \right \}  + j\left \{ e^{-j\frac{\pi}{4}n} + e^{j\frac{\pi}{4}n} \right \} \right ]u[n] $$
%
% $$ x[n]=\left | 0.8\right |^n  \left [ cos{\left \{  \frac{\pi n}{4}\right \}} + 2sin{\left \{\frac{\pi n}{4}\right \}} \right ]u[n] $$
%

%%
%verify in Matlab
n=0:20;
impseq=[1 zeros(1, length(n)-1)];         % generate impulse sequence
x1=filter(b,a,impseq);                    % compute impulse response
x2=(0.8.^n).*(cos(pi*n/4)+2*sin(pi*n/4)); % compute impulse response (book solution)

% Plot 
figure; subplot(211);
stem(n,x1); % plot impulse response
title('Impulse Response from filter fxn');
xlabel('n (samples)'); ylabel('Amplitude');

subplot(212);
stem(n,x2); % plot impusle response
title('Impulse Response from solution');
xlabel('n (samples)'); ylabel('Amplitude');

figure;
zplane(b,a); % plot PZ map


%% Problem 3: PFE and PZ-plot with repeated roots
% (Oppenheim, 3.33a)
%
% Find the partial fraction expansion of X(z) that is given by
%
% $$ X(z)=\frac {1} {(1-\frac{1}{2}z^{-1})^2(1-2z^{-1})(1-3z^{-1})} $$
% 
b=[1]; % numerator coefficients
a1=conv([1 1/2],[1 1/2]); % multiplying polynomials to get denom. coef.
a23=conv([1 -2],[1 -3]);
a=conv(a1,a23);

[r,p,k]=residuez(b,a) % compute residues (r), poles (p), and direct terms (k)

%%
% This represents the partial fraction expansion
%
% $$ X(z)=\frac{2.2041}{1-3z^{-1}}+\frac{-1.28}{1-2z^{-1}} + \frac{0.473}{1+\frac{1}{2}z^{-1}} + \frac{0.286}{(1+\frac{1}{2}z^{-1})^2} $$
%
figure; zplane(b,a)

%% Problem 4: Impulse, Step, and Freq. Response 
% (Oppenheim, Problem 2.9)
%
% Find the impulse response, freq. response, and step response of the
% LTI system defined by the difference equation
%
% $$ y[n]-\frac{5}{6}y[n-1]+\frac{1}{6}y[n-2]=\frac{1}{3}x[n-1]. $$
%

a=[1 -5/6 1/6];
b=[0 1/3 0];
figure; freqz(b,a) % frequency response
title('Freq. Response')

figure;
subplot(221); impz(b,a)  % plot impulse response
subplot(222); stepz(b,a) % plot step response

% Compare with book solution
n=0:14;
hImp=-2*(1/3).^n+2*(1/2).^n;
hStep=1+(1/3).^n-2*(1/2).^n;
subplot(223); stem(n,hImp)  % plot impulse response
title('Impulse Response (book)'); xlabel('n (samples)'); ylabel('Amplitude');
subplot(224); stem(n,hStep)  % plot impulse response
title('Step Response (book)'); xlabel('n (samples)'); ylabel('Amplitude');

%% Problem 5: Impulse of Causal and Anticausal System 
% (Oppenheim, Problem 2.16)
%
% Consider the difference equation
%
% $$ y[n]-\frac{1}{4}y[n-1]+\frac{1}{8}y[n-2]=3x[n-1] $$
%
% b) Both a causal and an anticausal LTI system are
% characterized by this difference equation. Find the impulse
% responses of the two systems.
% c) Show that the causal LTI is stable and the anticausal LTI
% system is unstable.
%
% From the difference equation and the ZT we can generate the system
% function H(z)
%
% $$ H(z)=\frac{3}{1-\frac{1}{4}z^{-1}-\frac{1}{4}z^{-1}} $$
%
% We can use matlab to get the expanded form of H(z).

b=[3];
a=[1 -1/4 -1/8];
[r,p,c]=residuez(b,a)

%%
% $$ H(z)=\frac{2}{1-0.5z^{-1}} + \frac{1}{1+0.25z^{-1}} $$
%
% Using Table 3.1 we get the following causal and anticausal
% impulse responses by taking the iZT of the terms in H(z).
%
% $$ h_c[n]=\left [ 2(1/2)^n + (-1/4)^n \right ] u[n] $$
% 
% $$ h_ac[n]=-\left [ 2(1/2)^n + (-1/4)^n \right ] u[-n-1] $$
%
% From ROC we see causal system is stable (ROC includes unit
% cirlce) and the anticausal system in unstable.

n=0:50;
hc=2*(1/2).^n + (-1/4).^n;
n=-50:0;
hac=-2*(1/2).^-n - (-1/4).^-n;
    
roc_c=max(abs(p));  %ROC of causal is righ-sided extending from largest pole 
roc_ac=min(abs(p)); %ROC of causal is righ-sided extending from smallest pole 

% Plot Pole-Zero maps with ROC
% Causal System
figure;
circle(roc_c,0) % draw ROC
hold('on');
zplane(b,a) % draw ZP map
title('PZ Map (causal)')
% Anticausal System
figure;
circle(roc_ac,0) % draw ROC
hold('on');
zplane(b,a) % draw ZP map
title('PZ Map (anticausal)')

%% Problem 5: Pole-Zero Plot
% (Oppenheim, Problem 3.12)
%
% Sketch the pole-zero plot for each of the following z-transforms and shade ROC:
%
% $$ X_1(z)=\frac{1-1/2z^{-1}}{1+1/2z^{-1}}, ROC: \left | z \right | <2 $$
%
% $$ X_2(z)=\frac{1-1/3z^{-1}}{\left ( 1+1/2z^{-1} \right )\left ( 1-2/3z^{-1} \right )}, x_2[n] causal $$
%
% $$ X_3(z)=\frac{1+z^{-1}-2z^{-2}}{1-13/6z^{-1}+z^{-2}} , x_3[n] absolutely summable $$
%

% define a and b of each z-function
a1=[1 2];
b1=[1 -1/2];
a2=conv([1 1/2],[1 -2/3]);
b2=[1 -1/3];
a3=[1 -13/6 1];
b3=[1 1 -2];

% find poles for each to determine ROC
p1=roots(a1);
p2=roots(a2);
p3=roots(a3);
disp(['p1=[' num2str(p1') '], p2=[' num2str(p2') '], p3=[' num2str(p3') ']' ])

%%
% X1(z) has ROC |z|<2, X2(z) has ROC |z|>2/3, X3(z) has ROC 2/3 < |z| < 3/2

%plot each PZ map and ROC
% X1
figure;
circle(p1,0) % draw ROC
hold 'on';
zplane(b1,a1) % draw ZP map
title('PZ Map: X_1(z)')

% X2
figure;
circle(max(abs(p2)),1) % draw ROC
hold 'on';
zplane(b2,a2); % draw ZP map
title('PZ Map: X_2(z)')
set(gca,'color',[.95 .95 .95])

% X3
figure;
circle(3/2,0); hold on; circle(2/3,2); % draw ROC
hold on;
zplane(b3,a3) % draw ZP map
title('PZ Map: X_3(z)')

%% References
% 
% Alan V. Oppenheim and Ronald W. Schafer. 2009. Discrete-Time
% Signal Processing (3rd ed.). Prentice Hall Press, Upper Saddle River, NJ, USA. 
% 
% Vinay K. Ingle and John G. Proakis. 1999. Digital Signal
% Processing Using MATLAB (1st ed.). Brooks/Cole Publishing Co., Pacific Grove, CA, USA.

%% Subfunctions

type circle
%%
type stepseq

%% Useful MATLAB Functions
%
% <html>
%  <ul>
%   <li><a href="http://www.mathworks.com/help/techdoc/ref/poly.html">poly</a></li>
%   <li><a href="http://www.mathworks.com/help/techdoc/ref/roots.html">roots</a></li>
%   <li><a href="http://www.mathworks.com/help/toolbox/signal/residuez.html">residuez</a></li>
%   <li><a href="http://www.mathworks.com/help/toolbox/control/ref/tf.html">tf</a></li>
%   <li><a href="http://www.mathworks.com/help/toolbox/control/ref/zpk.html">zpk</a></li>
%   <li><a href="http://www.mathworks.com/help/toolbox/signal/zp2tf.html">zp2tf</a></li>
%   <li><a href="http://www.mathworks.com/products/signal/functionlist.html">tf2zpk</a></li>
%   <li><a href="http://www.mathworks.com/help/techdoc/ref/conv.html">conv</a></li
%   <li><a href="http://www.mathworks.com/help/techdoc/ref/deconv.html">deconv</a></li>
%   <li><a href="http://www.engin.umich.edu/group/ctm/extras/commands.html">dimpulse</a></li>
%   <li><a href="http://www.mathworks.com/help/toolbox/signal/freqz.html">freqz</a></li>
%   <li><a href="http://www.mathworks.com/help/techdoc/ref/filter.html">filter</a></li>
%   <li><a href="http://www.mathworks.com/help/toolbox/control/ref/pzmap.html">pzmap</a></li>
%   <li><a href="http://www.mathworks.com/help/toolbox/signal/impz.html">impz</a></li>
%   <li><a href="http://www.mathworks.com/help/toolbox/signal/stepz.html">stepz</a></li>
%   <li><a
%   href="http://www.mathworks.com/help/techdoc/ref/poly.html">abc</a></li>
%  </ul>
% </html>

%%
close all;  % close all figures
clear       % clear workspace

