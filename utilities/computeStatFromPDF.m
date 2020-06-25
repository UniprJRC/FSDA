function out = computeStatFromPDF(Namefilepdf, sessione)
% FSDA undocumented to compute final mark from pdf input file

% second argument session can be E= estiva, or I= invernale or P =
% primaverile or M= magistrale

if nargin<2
    sessione='E';
end

% extract matricola and titolo tesi

% create a vote in numbers and in letters
votenum=[0 66:110]';
voteletters={'zero', 'sessantasei', 'sessantasette', 'sessantotto', 'sessantanove', ...
'settanta', 'settantuno', 'settantadue', 'settantatré', 'settantaquattro', ...
'settantacinque', 'settantasei', 'settantasette', 'settantotto', 'settantanove', ...
'ottanta', 'ottantuno', 'ottantadue', 'ottantatré', 'ottantaquattro', ...
'ottantacinque', 'ottantasei', 'ottantasette', 'ottantotto', 'ottantanove', ...
'novanta', 'novantuno', 'novantadue', 'novantatré', 'novantaquattro', ...
'novantacinque', 'novantasei', 'novantasette', 'novantotto', 'novantanove', ...
'cento', 'centouno', 'centodue', 'centotré', 'centoquattro', 'centocinque', ...
'centosei', 'centosette', 'centootto', 'centonove', 'centodieci'}';



str = extractFileText(Namefilepdf);

startmatr = strfind(str,"Matricola :");
endmatr = startmatr+6;

startnam = strfind(str,"Il Sig.");
endnam = strfind(str,"nato il");
startnaf = strfind(str,"La Sig.ra");
endnaf = strfind(str,"nata il");

startcorso = strfind(str, 'Ultimo piano di studio seguito ed anni accademici di iscrizion');
endcorso = strfind(str, "Insegnamento In Piano ");
endcorso=endcorso(1:2:length(endcorso)-1);

startrel = strfind(str, "Primo relatore:");
endrel = strfind(str, "Tipo della tesi:");
numtrenta = strfind(str, "30/30");
startlodi = strfind(str, "N. di lodi:");
endlodi = strfind(str, "N. Crediti Sovrannumerari");


startna=sort([startnam+8 startnaf+10]);
endna=sort([endnam endnaf]);

% totstr=sort([startna, numtrenta endna]);

starti = strfind(str,"Totale base 110");
endi = strfind(str,"Punti tesi:");

tabrows=numel(starti);
varNames = {'Matricola', 'Corso', 'Candidato','relatore','media base 110','Num di 30/30', 'Num lodi', 'InCorso' 'Voto_110'};
varTypes = {'string', 'string','string','string','double', 'double', 'double' 'string' 'double'};

n=numel(starti);
laureandi=table('Size',[tabrows numel(varTypes)], 'VariableTypes',varTypes, 'VariableNames',varNames);
jj=1;
for j=1:n
    
    % matricola
    start = startmatr(j);
    fin = endmatr(j);
    laureandi{j,1}=extractBetween(str,start-9,fin-20);
    
    % corso di laurea
    start = startcorso(j);
    fin = endcorso(j);
    AA=extractBetween(str,start,fin);
    %  findSlash=strfind(AA,"/");
    %    posLastSlash=findSlash(end);
    aa=AA{:};
    % Controlla se fuori corso
    if contains(aa,'FC')
        laureandi{j,7}="FC";
    else
        laureandi{j,7}="IC";
    end
    
    findCL=strfind(AA,"Corso di Laurea in");
    findCLsel= aa(findCL(end)+19:end-2);
    % findCR=regexp(findCLsel,'\n');
    findCLsel=removeExtraSpacesLF(findCLsel);
    laureandi{j,1}=string(findCLsel);
    
    % anagrafica candidato
    start = startna(j);
    fin = endna(j);
    laureandi{j,2}=extractBetween(str,start,fin-3);
    
    % relatore
    start = startrel(j);
    fin = endrel(j);
    laureandi{j,3}=extractBetween(str,start+16,fin-4);
    
    
    % media in 110
    start = starti(j);
    fin = endi(j);
    laureandi{j,4}=extractBetween(str,start+16,fin-2);
    
    % numero di trenta e 30/lode all'interno di un curriculum
    if  jj <= numel(numtrenta)
        trenta=0;
        
        while numtrenta(jj) < starti(j)
            trenta=trenta+1;
            jj=jj+1;
            if jj > numel(numtrenta)
                break
            end
        end
    end
    laureandi{j,5}=trenta;
    
    % numero di lodi
    start = startlodi(j);
    fin = endlodi(j);
    laureandi{j,6}=extractBetween(str,start+12,fin-4);
    
end

if strcmp(sessione,'E')
    perc=4;
elseif strcmp(sessione,'I')
    perc=3.5;
elseif strcmp(sessione,'P')
    perc=0;
elseif strcmp(sessione,'M')
    perc=0;
    magistrale=true;
else
    error('FSDA:computeStatFromPDF:WngIinput','Sessione deve essere E I oppure P')
end

% Premio percorso
AddInCorso=zeros(n,1);
AddInCorso(laureandi.InCorso=='IC')=perc;
if magistrale == true
    AddLodi=zeros(n,1);
    AddLodi(laureandi.("Num lodi")==2)=0.5;
    AddLodi(laureandi.("Num lodi")==3)=1;
    AddLodi(laureandi.("Num lodi")==4)=1;
    AddLodi(laureandi.("Num lodi")>=5)=1.5;    
else
    AddLodi=zeros(n,1);
    AddLodi(laureandi.("Num di 30/30")==2)=1;
    AddLodi(laureandi.("Num di 30/30")==3)=2;
    AddLodi(laureandi.("Num di 30/30")==4)=3;
    AddLodi(laureandi.("Num di 30/30")==5)=4;
    AddLodi(laureandi.("Num di 30/30")>5)=5;
end    
    laureandi.Voto_110=laureandi.("media base 110")+AddInCorso+AddLodi;


if magistrale==true
    % write table into Excel file
    filename = ['laureati-magistrali.xlsx'];
else
    filename = ['laureati-triennali.xlsx'];
end

writetable(laureandi,filename,'Sheet',1,'Range','A1')

out=laureandi;
end