clear all
set(0,'defaultfigurecolor','w')
% parpool(5)
tic

% nfp=1  stell-sym
% data_v=read_vmec('wout_fit3.nc');  
% data_b=read_boozer('boozmn_fit3.nc');

discrete_m=180;
discrete_n=180;
theta=linspace(0,2*pi,discrete_m+1); theta(end)=[]; theta=single(theta);
zeta=linspace(0,2*pi,discrete_n+1); zeta(end)=[]; zeta=single(zeta);


% filename='wout_fit3_ns=241_mpol=11_tol=-14.nc';
% filename='wout_2b40Ncoil16_B01opt_free.nc';
filename='wout_D3D.nc';
xm=ncread(filename,'xm')'; xm=single(xm);
xn=ncread(filename,'xn')'; xn=single(xn);
xm_nyq=ncread(filename,'xm_nyq')'; xm_nyq=single(xm_nyq);
xn_nyq=ncread(filename,'xn_nyq')'; xn_nyq=single(xn_nyq);
ns=ncread(filename,'ns');
iota=ncread(filename,'iotas')'; iota=reshape(iota,1,1,ns); iota=single(iota);
bsubumnc=ncread(filename,'bsubumnc'); bsubumnc=single(bsubumnc);
bsubvmnc=ncread(filename,'bsubvmnc'); bsubvmnc=single(bsubvmnc);
bsubsmns=ncread(filename,'bsubsmns'); bsubsmns=single(bsubsmns);
lmns=ncread(filename,'lmns'); lmns=single(lmns);
bmnc=ncread(filename,'bmnc'); bmnc=single(bmnc);
phi=ncread(filename,'phi')'; phi=single(phi);
chi=ncread(filename,'chi')'; chi=single(chi);


lambda=zeros(discrete_m,discrete_n,ns,'single');
d_lambda_d_theta=zeros(discrete_m,discrete_n,ns,'single');
d_lambda_d_zeta=zeros(discrete_m,discrete_n,ns,'single');
w=zeros(discrete_m,discrete_n,ns,'single');
d_w_d_theta=zeros(discrete_m,discrete_n,ns,'single');
d_w_d_zeta=zeros(discrete_m,discrete_n,ns,'single');
p=zeros(discrete_m,discrete_n,ns,'single');
d_p_d_theta=zeros(discrete_m,discrete_n,ns,'single');
d_p_d_zeta=zeros(discrete_m,discrete_n,ns,'single');
I=zeros(1,1,ns,'single');
G=zeros(1,1,ns,'single');

d_l_d_theta_mnc=zeros(size(lmns),'single');
d_l_d_zeta_mnc=zeros(size(lmns),'single');
wmns=zeros(size(bsubvmnc),'single');
d_w_d_theta_mnc=zeros(size(wmns),'single');
d_w_d_zeta_mnc=zeros(size(wmns),'single');
for i=1:size(lmns,1)
    d_l_d_theta_mnc(i,:)=lmns(i,:).*xm(1,i);
    d_l_d_zeta_mnc(i,:)=-lmns(i,:).*xn(1,i);
end
for i=1:size(bsubumnc,1)
    if (xm_nyq(1,i)~=0)
        wmns(i,:)=bsubumnc(i,:)./xm_nyq(1,i);
    elseif (xn_nyq(1,i)~=0)
        wmns(i,:)=-bsubvmnc(i,:)./xn_nyq(1,i);
    else 
        I(1,1,:)=bsubumnc(i,:);
        G(1,1,:)=bsubvmnc(i,:);
    end
    d_w_d_theta_mnc(i,:)=wmns(i,:).*xm_nyq(1,i);
    d_w_d_zeta_mnc(i,:)=-wmns(i,:).*xn_nyq(1,i);
end
clear bsubumnc bsubvmnc

mt=xm'*theta; 
mt_nyq=xm_nyq'*theta;
nz=xn'*zeta; 
nz_nyq=xn_nyq'*zeta; 
cosmt=cos(mt);  
cosmt_nyq=cos(mt_nyq);
sinmt=sin(mt);  
sinmt_nyq=sin(mt_nyq);
cosnz=cos(nz);  
cosnz_nyq=cos(nz_nyq);
sinnz=sin(nz);  
sinnz_nyq=sin(nz_nyq);
temp_xm=size(xm,2);
temp_xm_nyq=size(xm_nyq,2);
clear mt mt_nyq nz nz_nyq

for i=1:ns
    lmns_re=repmat(lmns(:,i),[1 discrete_m]);  
    d_l_d_theta_mnc_re=repmat(d_l_d_theta_mnc(:,i),[1 discrete_m]);  
    d_l_d_zeta_mnc_re=repmat(d_l_d_zeta_mnc(:,i),[1 discrete_m]);  
    wmns_re=repmat(wmns(:,i),[1 discrete_m]);  wmns_re=single(wmns_re);
    d_w_d_theta_mnc_re=repmat(d_w_d_theta_mnc(:,i),[1 discrete_m]);  
    d_w_d_zeta_mnc_re=repmat(d_w_d_zeta_mnc(:,i),[1 discrete_m]);  

    lambda(:,:,i)=(lmns_re.*sinmt)'*cosnz-(lmns_re.*cosmt)'*sinnz;
    d_lambda_d_theta(:,:,i)=(d_l_d_theta_mnc_re.*cosmt)'*cosnz+(d_l_d_theta_mnc_re.*sinmt)'*sinnz;
    d_lambda_d_zeta(:,:,i)=(d_l_d_zeta_mnc_re.*cosmt)'*cosnz+(d_l_d_zeta_mnc_re.*sinmt)'*sinnz;

    w(:,:,i)=(wmns_re.*sinmt_nyq)'*cosnz_nyq-(wmns_re.*cosmt_nyq)'*sinnz_nyq;  
    d_w_d_theta(:,:,i)=(d_w_d_theta_mnc_re.*cosmt_nyq)'*cosnz_nyq+(d_w_d_theta_mnc_re.*sinmt_nyq)'*sinnz_nyq;
    d_w_d_zeta(:,:,i)=(d_w_d_zeta_mnc_re.*cosmt_nyq)'*cosnz_nyq+(d_w_d_zeta_mnc_re.*sinmt_nyq)'*sinnz_nyq;
    
    p(:,:,i)=(w(:,:,i)-I(1,1,i).*lambda(:,:,i))./(G(1,1,i)+iota(1,1,i).*I(1,1,i));
    d_p_d_theta(:,:,i)=(d_w_d_theta(:,:,i)-I(1,1,i).*d_lambda_d_theta(:,:,i))./(G(1,1,i)+iota(1,1,i).*I(1,1,i));
    d_p_d_zeta(:,:,i)=(d_w_d_zeta(:,:,i)-I(1,1,i).*d_lambda_d_zeta(:,:,i))./(G(1,1,i)+iota(1,1,i).*I(1,1,i));
end
clear lmns d_l_d_theta_mnc d_l_d_zeta_mnc wmns d_w_d_theta_mnc d_w_d_zeta_mnc 
clear lmns_re d_l_d_theta_mnc_re d_l_d_zeta_mnc_re wmns_re d_w_d_theta_mnc_re d_w_d_zeta_mnc_re 
clear cosmt sinmt cosnz sinnz d_w_d_theta d_w_d_zeta

theta_b=zeros(discrete_m,discrete_n,ns,'single');
zeta_b=zeros(discrete_m,discrete_n,ns,'single');
d_boozer_d_vemc=zeros(discrete_m,discrete_n,ns,'single');
for i=1:discrete_m
    for j=1:discrete_n
        theta_b(i,j,:)=theta(1,i)+lambda(i,j,:)+iota(1,1,:).*p(i,j,:);
        zeta_b(i,j,:)=zeta(1,j)+p(i,j,:);
        d_boozer_d_vemc(i,j,:)=(1+d_lambda_d_theta(i,j,:)).*(1+d_p_d_zeta(i,j,:))+d_p_d_theta(i,j,:).*(iota(1,1,:)-d_lambda_d_zeta(i,j,:));
    end
end
clear lambda p d_lambda_d_theta d_lambda_d_zeta d_p_d_theta d_p_d_zeta w
B_v=zeros(discrete_m,discrete_n,ns,'single');
Bs_v=zeros(discrete_m,discrete_n,ns,'single');
for i=1:ns
    bmnc_re=repmat(bmnc(:,i),[1 discrete_m]);  bmnc_re=single(bmnc_re);
    bsubsmns_re=repmat(bsubsmns(:,i),[1 discrete_m]);  bsubsmns_re=single(bsubsmns_re);
    B_v(:,:,i)=(bmnc_re.*cosmt_nyq)'*cosnz_nyq+(bmnc_re.*sinmt_nyq)'*sinnz_nyq;
    Bs_v(:,:,i)=(bsubsmns_re.*sinmt_nyq)'*cosnz_nyq-(bsubsmns_re.*cosmt_nyq)'*sinnz_nyq; 
end
clear bmnc bsubsmns bmnc_re bsubsmns_re cosmt_nyq sinmt_nyq cosnz_nyq sinnz_nyq
toc

mpol=max(xm);
ntor=max(xn);
sizem=mpol+1;
sizen=2*ntor+1;
B_mnc_b=zeros(sizem,sizen,ns,'single');
Bs_mns_b=zeros(sizem,sizen,ns,'single');
delta=(theta(1,2)-theta(1,1))*(zeta(1,2)-zeta(1,1));
temp_B_v=delta*d_boozer_d_vemc.*B_v; temp_B_v=single(temp_B_v);
temp_Bs_v=delta*d_boozer_d_vemc.*Bs_v; temp_Bs_v=single(temp_Bs_v);
clear Bs_v B_v d_boozer_d_vemc

B_mnc_b_re=zeros(size(xm,2),ns,'single');
Bs_mns_b_re=zeros(size(xm,2),ns,'single');
cosm=zeros(max(xm)+1,size(theta_b,1),size(theta_b,2),size(theta_b,3),'single');
sinm=zeros(max(xm)+1,size(theta_b,1),size(theta_b,2),size(theta_b,3),'single');
cosn=zeros(max(xn)+1,size(zeta_b,1),size(zeta_b,2),size(zeta_b,3),'single');
sinn=zeros(max(xn)+1,size(zeta_b,1),size(zeta_b,2),size(zeta_b,3),'single');
    
for j=1:temp_xm   % parfor

    if (xm(1,j)==0 && xn(1,j)<0)
        B_mnc_b_re(j,:)=0;
        Bs_mns_b_re(j,:)=0;
    else
        B_mnc_b_re(j,:)=sum(temp_B_v.*cos(xm(1,j)*theta_b-xn(1,j)*zeta_b),[1 2])/(2*pi*pi);
        Bs_mns_b_re(j,:)=sum(temp_Bs_v.*sin(xm(1,j)*theta_b-xn(1,j)*zeta_b),[1 2])/(2*pi*pi);           
    end
end
clear temp_B_v temp_Bs_v theta_b zeta_b

for j=1:temp_xm
    m=xm(j)+1;
    n=xn(j)+ntor+1;
    B_mnc_b(m,n,:)=B_mnc_b_re(j,:);
    Bs_mns_b(m,n,:)=Bs_mns_b_re(j,:);
end
B_mnc_b(1,ntor+1,:)=B_mnc_b(1,ntor+1,:)/2;
clear B_mnc_b_re Bs_mns_b_re 
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear i j k m n x y bcedge aa aaa k1 k2 selmp selnt bmnmax b00 bmn nfs rho rowmpol rowntor apbc
rowmpol=0:mpol;                     
rowntor=-1*ntor:ntor;    

bmn=B_mnc_b;

nfs=ns;
rho=ncread(filename,'phi')';
rho=(rho/rho(length(rho))).^0.5;
rho=rho(1:nfs);



% bcedge=abs(bmn(:,:,21));
 for temp1=1:size(bmn,1)
     for temp2=1:size(bmn,2)
         for temp3=1:size(bmn,3)-(ns-1)/10*7
          %   bcedge(temp1,temp2)=abs(bmn(temp1,temp2,temp3))*(rho(temp3)-rho(temp3-1));
            bcedge(temp1,temp2)=abs(bmn(temp1,temp2,temp3))*rho(temp3);
         end
     end
 end
 clear temp1 temp2 temp3



 [x,y]=find(bcedge==max(max(bcedge)));
 b00=squeeze(bmn(x,y,:))';
 bcedge(x,y)=0;
 
 
 i=0;
 figure
for k=1:10
apbc=max(max(bcedge));
[x,y]=find(bcedge==apbc);
x=x(1);
y=y(1);
i=i+1;
m(i)=x;
n(i)=y;
bcedge(x,y)=0;
selmp(i)=rowmpol(x);
selnt(i)=rowntor(y);
bmnmax(i,:)=squeeze(bmn(x,y,:));
k1=bmnmax(i,1:nfs);
k2=k1(2:nfs)./b00(2:nfs);
% k2=k1(2:nfs);
k2=[0 k2];
plot(rho,k2,'LineWidth',4);
hold on
end

aa=arrayfun(@(x,y)sprintf('%.0f,',x,y),selmp,selnt,'uniformoutput',0);
aaa=char(aa);
for temp1=1:1:size(aaa,1)
    for temp2=1:1:10
    if aaa(temp1,size(aaa(temp1,:),2)-temp2+1)==(',')
       aaa(temp1,size(aaa(temp1,:),2)-temp2+1)=(' ');
        break;
    end
    end
    clear temp2
end
clear temp1

    legend(aaa,'FontSize',25,'FontWeight','bold')
%     xlabel('\rho_ ','FontSize',30,'FontWeight','bold'),ylabel('B_m_,_n','FontSize',30,'FontWeight','bold')
    xlabel('\rho_ ','FontSize',30,'FontWeight','bold'),ylabel('B_m_,_n/B_0_,_0','FontSize',30,'FontWeight','bold')
    xlim([0 rho(nfs)]);
%     ylim([-1 1]);
%     ylim([-0.2 0.2]);
    xticks([0 0.25 0.5 0.75 1])
     set(gca,'FontSize',40,'FontWeight','bold','linewidth',3);
     box on;
     axis square
     clear i j k m n x y bcedge aa aaa k1 k2 selmp selnt bmnmax b00 bmn nfs rho rowmpol rowntor apbc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear i j k m n x y bcedge aa aaa k1 k2 selmp selnt bmnmax b00 bmn nfs rho rowmpol rowntor apbc
rowmpol=0:mpol;                     
rowntor=-1*ntor:ntor;    

bmn=Bs_mns_b;

nfs=ns;
rho=ncread(filename,'phi')';
rho=(rho/rho(length(rho))).^0.5;
rho=rho(1:nfs);

bcedge=abs(bmn(:,:,21));
 for temp1=1:size(bmn,1)
     for temp2=1:size(bmn,2)
         for temp3=1:size(bmn,3)-(ns-1)/10*7
          %   bcedge(temp1,temp2)=abs(bmn(temp1,temp2,temp3))*(rho(temp3)-rho(temp3-1));
            bcedge(temp1,temp2)=abs(bmn(temp1,temp2,temp3))*rho(temp3);
         end
     end
 end
 clear temp1 temp2 temp3

 
 i=0;
 figure
for k=1:10
apbc=max(max(bcedge));
[x,y]=find(bcedge==apbc);
x=x(1);
y=y(1);
i=i+1;
m(i)=x;
n(i)=y;
bcedge(x,y)=0;
selmp(i)=rowmpol(x);
selnt(i)=rowntor(y);
bmnmax(i,:)=squeeze(bmn(x,y,:));
k1=bmnmax(i,1:nfs);
k2=k1(2:nfs);
% k2=k1(2:nfs);
k2=[0 k2];
plot(rho,k2,'LineWidth',4);
hold on
end

aa=arrayfun(@(x,y)sprintf('%.0f,',x,y),selmp,selnt,'uniformoutput',0);
aaa=char(aa);
for temp1=1:1:size(aaa,1)
    for temp2=1:1:10
    if aaa(temp1,size(aaa(temp1,:),2)-temp2+1)==(',')
       aaa(temp1,size(aaa(temp1,:),2)-temp2+1)=(' ');
        break;
    end
    end
    clear temp2
end
clear temp1

    legend(aaa,'FontSize',25,'FontWeight','bold')
    xlabel('\rho_ ','FontSize',30,'FontWeight','bold'),ylabel('B_{S} _m_,_n','FontSize',30,'FontWeight','bold')
%     xlabel('\rho_ ','FontSize',30,'FontWeight','bold'),ylabel('B_m_,_n/B_0_,_0','FontSize',30,'FontWeight','bold')
    xlim([0 rho(nfs)]);
%     ylim([-1 1]);
%     ylim([-0.2 0.2]);
    xticks([0 0.25 0.5 0.75 1])
     set(gca,'FontSize',40,'FontWeight','bold','linewidth',3);
     box on;
     axis square
clear i j k m n x y bcedge aa aaa k1 k2 selmp selnt bmnmax b00 bmn nfs rho rowmpol rowntor apbc


