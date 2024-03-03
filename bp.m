radky=100;
iterace = 100; %počet sloupců = počet řádků + i, i=1,...,iterace
hustota = 0.1; %hustota řídké matice
osa_x = 1:iterace;

chyby_fr = zeros(4,iterace); %chyba a čas obyčejné náhodné matice
casy_fr = zeros(4,iterace);

chyby_sr = zeros(4,iterace);%chyba a čas řídké matice
casy_sr = zeros(4,iterace);

for i=radky+1:radky+iterace %zpracování dat - vygenerujeme pro každou iteraci nové matice, aplikujeme na něj všechny 4 metody a uložíme naměřené chyby a časy výpočtů do příslušných matic dat
    fullrandom = rand(radky,i); %full náhodná matice
    sparserandom = sprand(radky,i,hustota); %řídká náhodná matice
    m=radky;
    n=i;
    [chyba_fr,cas_fr] = KernelCalculation(fullrandom,m,n);
    [chyba_sr,cas_sr] = KernelCalculation(sparserandom,m,n);

    chyby_fr(:,i-radky) = chyba_fr;
    casy_fr(:,i-radky) = cas_fr;
    
    chyby_sr(:,i-radky) = chyba_sr;
    casy_sr(:,i-radky) = cas_sr;  
end
figure(1);
p1=plotdata(osa_x,chyby_fr,'Norma chyby výpočtu pro náhodnou plnou matici','počet sloupců - počet řádků','Norma chyby');
figure(2);
p2=plotdata(osa_x,casy_fr,'Čas výpočtu pro náhodnou plnou matici','počet sloupců - počet řádků','Čas výpočtu');
figure(3);
p3=plotdata(osa_x,chyby_sr,'Norma chyby výpočtu pro náhodnou řídkou matici','počet sloupců - počet řádků','Čas výpočtu');
figure(4);
p4=plotdata(osa_x,casy_sr,'Čas výpočtu pro náhodnou řídkou matici','počet sloupců - počet řádků','Čas výpočtu');

function [norma,time] = LUfactorization(A,m,n) %Výpočet LU
    tic

    [~,U] = lu(A);
    
    U1=U(:,1:m);
    U2=U(:,m+1:n);
    J=inv(U1)*U2;
    B=[-J;eye(n-m)];

    norma = norm(A*B,"fro");
    time=toc;
end

function [norma,time] = GJfactorization(A,m,n) %výpočet Gauss-Jordana
    tic

    R=rref(A);
    J=R(:,m+1:n);
    B=[-J;eye(n-m)];

    norma=norm(A*B,"fro");
    time=toc;
end

function [norma,time] = QRfactorization(A,m,n) %výpočet QR
    tic

    [~,R]=qr(A);
    R1=R(:,1:m);
    R2=R(:,m+1:n);
    J=inv(R1)*R2;
    B=[-J;eye(n-m)];
    
    norma=norm(A*B,"fro");
    time=toc;
end

function [norma,time]= SVDfactorization(A,m,n) %výpočet SVD
    tic

    [~,~,V] = svd(full(A));
    Q=V';
    Q2=Q(:,m+1:n);
    B=Q2;

    norma=norm(A*B,"fro");
    time=toc;
end

function [normvector,timevector] = KernelCalculation(A,m,n) %Vrací chybový a časový vektor výpočtů všech metod pro jednu stejnou matici A
    [n1,t1] = LUfactorization(A,m,n);
    [n2,t2] = GJfactorization(A,m,n);
    [n3,t3] = QRfactorization(A,m,n);
    [n4,t4] = SVDfactorization(A,m,n);

    normvector = [n1 n2 n3 n4];
    timevector = [t1,t2,t3,t4];

end

function p = plotdata(x,y,title_text,x_axis_text,y_axis_text)
    p=plot(x,y);
    title(title_text)
    xlabel(x_axis_text)
    ylabel(y_axis_text,'Rotation',0)
    legend('LU','GJE','QR','SVD')
end