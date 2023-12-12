function [V,DD]= GPR_active_sub(myGPR)

% Construct C matrix
nmcs = 1e3;
nvar = size(myGPR.ExpDesign.X,2);
XRAND = uq_getSample(nmcs)*0.999;


grd = zeros(size(XRAND,1),nvar);
dd = (0.001);
for jj = 1:size(XRAND,1)
    x = XRAND(jj,:);
    xl = repmat(x,nvar,1);
    xu = repmat(x,nvar,1);
    for kk = 1:nvar
        xl(kk,kk) = x(kk)-dd;
        xu(kk,kk) = x(kk)+dd;
    end
    resu = uq_evalModel(myGPR,[xu]);
    resl = uq_evalModel(myGPR,[xl]);
    grd(jj,:) =(resu-resl)'./(2*dd)   ;
    jj
end

C = (grd'*grd)/nmcs;
[V D] = eig(C);
DD = flipud(diag(D));
V = fliplr(V);
DFG = DD/sum(DD);
cums = cumsum(DFG);
