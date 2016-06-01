function [A, X, Ag0, Ap, F0] = VCV_inputs(folder)
%takes folder as input (i.e., VCV/real2d/asa-area)

M = 100;

%open files
fileAF = strcat(folder, '/file.AF');
fileAg0 = strcat(folder, '/file.AG0');
fileAgp = strcat(folder, '/file.AGP');
fileF0 = strcat(folder, '/file.F0');
AFID = fopen(fileAF);
Ag0ID = fopen(fileAg0);
AgpID = fopen(fileAgp);
F0ID = fopen(fileF0);

%initialize arrays

Ag0data = textscan(Ag0ID, '%f %f %s');
length = Ag0data{1}(1);
[x, v] = deal(zeros(length, 1));
maxlength = Ag0data{1}(length+1)*10;

[Ag0, Ap, F0] = deal(zeros(maxlength, 1));
[A, X] = deal(zeros(maxlength, M));

%create Ag0 array
for i = 1:length
    x(i) = (Ag0data{1}(i+1)) * 10 + 1;
    v(i) = Ag0data{2}(i+1);
end

for i = 1:length-1
    s = Ag0data{3}(i+2);
    if strcmp(s, 'COS') == 1
        t = 'pchip';
    else t = 'linear';
    end
    Ag0(x(i):x(i+1)-1) = interp1(x, v, x(i):x(i+1)-1, t);
end

%create Agp array

Agpdata = textscan(AgpID, '%f %f %s');
length = Agpdata{1}(1);
[x, v] = deal(zeros(length, 1));

for i = 1:length
    x(i) = (Agpdata{1}(i+1)) * 10 + 1;
    v(i) = Agpdata{2}(i+1);
end

for i = 1:length-1
    s = Agpdata{3}(i+2);
    if strcmp(s, 'COS') == 1
        t = 'pchip';
    else t = 'linear';
    end
    Ap(x(i):x(i+1)-1) = interp1(x, v, x(i):x(i+1)-1, t);
end

%create F0 array

F0data = textscan(F0ID, '%f %f %s');
length = F0data{1}(1);
[x, v] = deal(zeros(length, 1));

for i = 1:length
    x(i) = (F0data{1}(i+1)) * 10 + 1;
    v(i) = F0data{2}(i+1);
end

for i = 1:length-1
    s = F0data{3}(i+2);
    if strcmp(s, 'COS') == 1
        t = 'pchip';
    else t = 'linear';
    end
    F0(x(i):x(i+1)-1) = interp1(x, v, x(i):x(i+1)-1, t);
end
F0(3000:min(5000, end)) = 100;

%create area function

width = str2double(fgetl(AFID));
length = str2double(fgetl(AFID));
Adata = zeros(width+1, length);
Xdata = zeros(width+1, length);
[Ai, Xid] = deal(zeros(width+1, 100));
x = transpose(1:100:(width+1)*100);
for i = 0:width
    areafile = strcat(folder, '/', num2str(i*10), '.ARE');
    areaID = fopen(areafile);
    title = str2double(fgetl(areaID));
    length = str2double(fgetl(areaID));
    xval = str2double(fgetl(areaID));
    Xdata(i+1,:) = xval;
    C = textscan(areaID, '%f');
    Adata(i+1,:) = C{1};
    
    Xc = cumsum(Xdata(i+1,:));
    Xi = linspace(Xc(1), Xc(end), M);
    Ai(i+1,:) = interp1(Xc, Adata(i+1,:), Xi);
    Xid(i+1,:) = diff([0 Xi]);
    
    fclose(areaID);
end

for j = 1:M
    Av = Ai(:,j);
    Xv = Xid(:,j);
    for i = 1:width-1
        A(x(i):x(i+1)-1,j) = interp1(x, Av, x(i):x(i+1)-1, 'pchip');
        X(x(i):x(i+1)-1,j) = interp1(x, Xv, x(i):x(i+1)-1, 'pchip');
    end
end

%close files
fclose(AFID);
fclose(Ag0ID);
fclose(AgpID);
fclose(F0ID);

end