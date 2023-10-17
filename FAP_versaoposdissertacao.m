%   Rotina de pre-processamento de sinais de ERP com base em filtragem adaptativa.

versao=version;
if IsOctave,
  pkg load signal;
  pkg load statistics;
end;
if versao(1)>'4', warning('off'); end;
BW=0; %Black & White plots
M=240;
Fs=400;
N=Fs;
NCh=4; %Fp1, Fp2, Cz e Oz
%Matriz de mistura (Fp1, Fp2, Cz e Oz):
Mix=sqrt([0.5 0.3 0.15 0.05; ...
    0.3 0.5 0.15 0.05; ...
    0.2 0.2 0.4  0.2;  ...
    0.2 0.2 0.2  0.4]);
global M
global B;
A=20; %fator de escala do ERP
A2=200; %fator de escala da piscada
Sigma=50; %d.p. do ru�do branco que gerar� o EEG de fundo
tpe=[0:N-1]/Fs;
B=5; Padrao=A*SimFunc(tpe); if size(Padrao,2)==1, Padrao=Padrao'; end;
Zn=[0.5;0.5;1;0.8]*Padrao;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coeficientes do filtro do EEG (Fp1, Fp2, Cz e Oz):
freqsup=[35 35 25 12];
freqinf=[ 8  8  5 0.5];
for i=1:NCh,
  [bleeg(i,:),aleeg(i,:)]=butter(2,2*freqsup(i)/Fs);
  [bheeg(i,:),aheeg(i,:)]=butter(2,2*freqinf(i)/Fs,'high');
  beeg(i,:)=conv(bheeg(i,:),bleeg(i,:));
  aeeg(i,:)=conv(aheeg(i,:),aleeg(i,:));
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Percentuais de �pocas em que v�o ocorrer piscadas
Percs=100;%20:30:80; %0:10:100;
%Par�metros da filtragem adaptativa:
LL=[21]; %ordens da filtragem adaptativa
%alpha=[1e-5];% 1e-6 1e-5]; %fator de adapta��o
alpha=[1e-7]
%alpha=[1e-7 1e-6 1e-5]; %fator de adapta��o
%Par�metros do peso para a filtragem adaptativa
RMSmins=[10]; %[10 30]; %30;
RMSmaxs=[60]; %[40 60]; %60;
[bk,ak]=butter(2,3/Fs*2); %passa-baixas em 1Hz para o RMS
%Tempos de piscada:
Tps=[0.25 0.5 0.75]; %0.25; %tempos m�dios de ocorr�ncia das piscadas;
Tv=0.25; %desvio-padr�o dos momentos das piscadas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Primeiro, gera o EEG simulado:
disp('Gerando o sinal EEG cru...');
Rn=randn(NCh,N*M)*Sigma;
for j=1:NCh,
  Rn(j,:)=filter(beeg(j,:),aeeg(j,:),Rn(j,:));
end;
EEGERP=Mix*Rn;
%Depois, soma os ERPs
for i=1:M,
  EEGERP(:,(i-1)*N+1:i*N)=Zn+EEGERP(:,(i-1)*N+1:i*N);
end;
EEG=zeros(size(EEGERP));
clear Rn;
%Analisa entre as v�rias ocorr�ncias das piscadas:
for iP=1:length(Percs),
  Perc=Percs(iP);
  ss=rand(1,M); ll=prctile(ss,Perc); Pisca=ss<=ll;
  for iT=1:length(Tps),
    Tp=Tps(iT);
    delay=(rand(1,M)-0.5)*Tv+Tp;
    for i=1:M,
      B=1; if Pisca(i), Piscada=SimBlink(tpe,delay(i)); else, Piscada=zeros(1,N); end;
      if size(Piscada,2)==1, Piscada=Piscada'; end;
      %%%%%%%%%%%%%%
      [bk_teste1,ak_teste1]=butter(2,7/Fs*2); %passa-baixas em 7Hz para o Cz
      [bk_teste2,ak_teste2]=butter(2,5/Fs*2); %passa-baixas em 5Hz para o Oz
      Piscadafilt_Cz=filter(bk_teste1,ak_teste1,Piscada);
      Piscadafilt_Cz= Piscadafilt_Cz * 1/2;
      %piscada em Oz apenas com filtro
      Piscadafilt_Oz=filter(bk_teste2,ak_teste2,Piscada);
      %piscada em Oz apenas com filtro e ponderada
      Piscadafilt_Oz = Piscadafilt_Oz * 1/4;     
      PiscadaOriginal_Oz = Piscada * 1/4;
      paramA = 0.8;
      paramB = 0.6;
      %paramB = 0.4;
      %piscada em Oz apenas não linear
      Piscada_naolinear_Oz = paramA*Piscada + paramB*(Piscada.^3);
      %piscada em Oz apenas não linear e ponderada
      Piscada_naolinear_Oz_ponderada = Piscada_naolinear_Oz * 1/4;
      %piscada em Oz com filtro e nao linear
      Piscadafilt_naolinear_Oz=filter(bk_teste2,ak_teste2,Piscada_naolinear_Oz);
      % %piscada em Oz com filtro e nao linear e ponderada
      Piscadafilt_naolinear_Oz = Piscadafilt_naolinear_Oz * 1/4;
      PiscadaTotal = [Piscada; Piscada; Piscadafilt_Cz; Piscadafilt_naolinear_Oz];
      Yn=A2*PiscadaTotal;
      %Yn=A2*[1;1;1/2;1/4]*Piscada;
      EEG(:,(i-1)*N+1:i*N)=Yn+EEGERP(:,(i-1)*N+1:i*N);
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RMS=(EEG(1,:)+EEG(2,:))/2;
    RMS=filtfilt(bk,ak,RMS.^2); RMS=sqrt(RMS);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iMin=1:length(RMSmins),
      RMSmin=RMSmins(iMin);
      for iMax=1:length(RMSmaxs),
        RMSmax=RMSmaxs(iMax);
        %Calcula o peso Kapa
        Kapa=(RMS-RMSmin)/(RMSmax-RMSmin);
        iLower=find(Kapa<0); iHigher=find(Kapa>1);
        Kapa(iLower)=zeros(size(iLower)); Kapa(iHigher)=ones(size(iHigher));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for iL=1:length(LL),
          L=LL(iL);
          for iAlfa=1:length(alpha),
            ERPa=zeros(size(Zn)); ERPb=ERPa; ERPc=ERPa;
            alfa=alpha(iAlfa);
            %Inicializacoes:
            FAdap=zeros(L,NCh);
            BuffAdap=zeros(NCh,L-1);
            if versao(1)>'4', h=waitbar(0,'Processando...'); else, disp(' '); end;
            CondIniEEG=zeros(max(size(aeeg,2),size(beeg,2))-1,NCh);
            SampAcum=0;
            for i=1:M,
              Xn=EEG(:,(i-1)*N+1:i*N);
              ERPa=ERPa+Xn;
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%
              Sinal=[BuffAdap,Xn];
              Ruido=(Sinal(1,:)+Sinal(2,:))/2;
              %Filtragem adaptativa propriamente dita:
              for j=1:N,
                X=Ruido(j+L-1:-1:j); %1 x L
                E=Sinal(:,j+L-1)'-X*FAdap; % 1 x NCh
                G=Sinal(:,j+L-1)'-X*FAdap*Kapa(SampAcum+j); % 1 x NCh
                FAdap=FAdap+alfa*(E'*X)'; %L x NCh
                Xn(:,j)=G';
              end;
              BuffAdap=Sinal(:,size(Sinal,2)-L+2:size(Sinal,2));
              SampAcum=SampAcum+N;
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%
              ERPb=ERPb+Xn;
              if versao(1)>'4',
                waitbar(i/M);
              else,
                fprintf(1,'.');
              end;
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%
              Wn=EEGERP(:,(i-1)*N+1:i*N);
              ERPc=ERPc+Wn;
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%
              %ERPa é o sinal com piscada sem filtro, o sinal contaminado
              %ERPb é o sinal após o procedimento de filtragem
              %ERPc é o valor verdadeiro (vetor de interesse, EEG com ERP antes de adicionar interferencia da piscada)
              % devemos comparar o ERPb com o ERPc
              
            end;
            if versao(1)>'4', close(h); else, disp(' '); end;
            ERPa=ERPa/M; ERPb=ERPb/M; ERPc=ERPc/M;
            iplot=[1,2,3,5];
            ChNames=['Fp1';'Fp2';'Cz ';'Oz '];
            figure; whitebg(gcf,'w'); set(gcf,'color','w','name',['L=',int2str(L),' - iAlfa=',int2str(iAlfa),' - Tp=',int2str(round(1000*Tp)),...
                ' - RMSmin=',int2str(RMSmin),' - RMSmax=',int2str(RMSmax),' - Perc=',int2str(Perc)]);
            for i=1:NCh,
              ax(i)=subplot(3,2,iplot(i));
              
                ll=plot(1000*tpe,ERPc(i,:),'g',1000*tpe,ERPa(i,:),'b-',1000*tpe,ERPb(i,:),'k-');
                set(ll(1),'linewidth',2);

              ylabel([deblank(ChNames(i,:)),' (',setstr(181),'V)']);
              set(gca,'ylim',[-1 1]*2*A);
              if versao(1)>'4',
                ErroRMS(iL,iAlfa,i,iT,iMin,iMax,iP)=sqrt(mean((ERPb(i,:)-ERPc(i,:)).^2))/sqrt(mean(ERPc(i,:).^2));
              else,
                eval(['ErroRMS_Ch',int2str(i),'_Tp',int2str(iT),'(iL,iAlfa)=sqrt(mean((ERPb(i,:)-ERPc(i,:)).^2));']);
              end;
            end;
            xlabel('Time (ms)');
            set(ax,'unit','normal');
            wd=get(ax(3),'pos'); wd(1)=(1-wd(3))/2; set(ax(3),'pos',wd);
            wd=get(ax(4),'pos'); wd(1)=(1-wd(3))/2; set(ax(4),'pos',wd);
            set(ax(1:3),'xticklabel','');
            if i==NCh,
              if BW,
                title(' ');
              else,
                title(' ');
                disp(Tp);
                disp(alfa);
                disp(L);
                MSE1 = immse(ERPb(1,:),ERPc(1,:));
                RMSE1 = sqrt(MSE1);
                disp(RMSE1);
                MSE2 = immse(ERPb(2,:),ERPc(2,:));
                RMSE2 = sqrt(MSE2);
                disp(RMSE2);
                MSE3 = immse(ERPb(3,:),ERPc(3,:));
                RMSE3 = sqrt(MSE3);
                disp(RMSE3);
                MSE4 = immse(ERPb(4,:),ERPc(4,:));
                RMSE4 = sqrt(MSE4);
                disp(RMSE4);
                MSEmedio = (MSE1+MSE2+MSE3+MSE4)/4;
                RMSEmedio = (RMSE1+RMSE2+RMSE3+RMSE4)/4;
                disp('valor medio');
                disp(RMSEmedio);
              end;
            end;        
          end; %for iAlfa
        end; % for iL LL length
      end; % for iMax RMSmax
    end; % for iMin RMSmin
  end; % for iT Tps
end;
if versao(1)>'4', save ResultAF ErroRMS; else, save ResultAF ErroRMS1 ErroRMS2 ErroRMS3 ErroRMS4; end;
