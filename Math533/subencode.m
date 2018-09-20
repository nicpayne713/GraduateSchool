function cipher=subencode(message)
capabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
alphabet='abcdefghijklmnopqrstuvwxyz';
for i=1:26,
    message=strrep(message,capabet(i),alphabet(i));
end
X=1;
while X<=length(message),
    if (double(message(X))<97)&&not(double(message(X))==32);
        message(X)='';
        X=X+1;
    elseif (double(message(X))>122)&&not(double(message(X))==32);
        message(X)='';
        X=X+1;
    else X=X+1;
    end
end
cipherbet=char(randperm(26)+64);
for i=1:26,
    message=strrep(message,alphabet(i),cipherbet(i));
end
for i=1:26,
cipherbet=strrep(cipherbet,capabet(i),alphabet(i));
message=strrep(message,capabet(i),alphabet(i));
end
disp('here is the encoded alphabet and message')
alphabet=cat(1,alphabet,cipherbet)
cipher=message;
