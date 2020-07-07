input = fopen('auxMm.txt', 'r');
outpt = fopen('memor.txt', 'w');
dt = fscanf(input, '%f', 1);
ofDp = fscanf(input, '%f', 1);
szDp = fscanf(input, '%f', 1);
pN = fscanf(input, '%i', 1);
pAsp = fscanf(input, '%f', pN);
pAss = fscanf(input, '%f', [pN,pN]);
pAss = transpose(pAss);
pN = pN + 1;
outpt = fopen('memor.txt', 'w');
for i = 1:szDp
    r = (ofDp / szDp) * (i - 0.5);
    lambda = sqrt(10500*exp(-7.38*r));
    pApp = zeros(pN,pN);
    pApp(1,2:pN) = -lambda * pAsp(:);
    pApp(2:pN,1) = lambda * pAsp(:);
    pApp(2:pN,2:pN) = pAss(:,:);
    pT = expm(-dt * pApp);
    pS = eye(pN) - pT * transpose(pT);
    pS = sqrtm(pS);
    1 / (pApp(1,1) - pApp(1,2:pN) * inv(pApp(2:pN,2:pN)) * pApp(2:pN,1))
    r
    pApp
    for i = 1:pN
        for j = 1:pN
            fprintf(outpt, '%2.8f\t ', pT(i,j));
        end
        fprintf(outpt, '\n');
    end
    fprintf(outpt, '\n');
    for i = 1:pN
        for j = 1:pN
            fprintf(outpt, '%2.8f\t ', pS(i,j));
        end
        fprintf(outpt, '\n');
    end
    fprintf(outpt, '\n');
end

oN = fscanf(input, '%i', 1);
oAsp = fscanf(input, '%f', oN);
oAss = fscanf(input, '%f', [oN,oN]);
oAss = transpose(oAss);
oN = oN + 1;
for i = 1:szDp
    r = (ofDp / szDp) * (i - 0.5);
    lambda = sqrt(12591*exp(-7.36*r));
    oApp = zeros(oN,oN);
    oApp(1,2:oN) = -lambda * oAsp(:);
    oApp(2:oN,1) = lambda * oAsp(:);
    oApp(2:oN,2:oN) = oAss(:,:);
    oT = expm(-dt * oApp);
    oS = eye(oN) - oT * transpose(oT);
    oS = sqrtm(oS);
    r
    oApp
    for i = 1:oN
        for j = 1:oN
            fprintf(outpt, '%2.8f\t ', oT(i,j));
        end
        fprintf(outpt, '\n');
    end
    fprintf(outpt, '\n');
    for i = 1:oN
        for j = 1:oN
            fprintf(outpt, '%2.8f\t ', oS(i,j));
        end
        fprintf(outpt, '\n');
    end
    fprintf(outpt, '\n');
end

%{
pT = zeros(pN,pN);
pTT = zeros(pN,pN);
pT(2:pN,2:pN) = expm(-dt * pApp(2:pN,2:pN));
pTT(2:pN,2:pN) = expm(-dt * transpose(pApp(2:pN,2:pN)));
pS = eye(pN) - pT * transpose(pT);
pS = sqrtm(pS);
pS(1,1) = 0;
dt
pT
pS
%}