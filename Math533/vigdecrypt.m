function attempt=vigdecrypt(message)
%This function estimates the key length for a message encrypted with a
%vigenere ciphter, then works to decode the message.
%It removes all punctuation and case from the original
%text.
tic
alphabet='abcdefghijklmnopqrstuvwxyz';
%Alphabet = ['A' , 'B' , 'C' , 'D' , 'E' , 'F' , 'G' , 'H', 'I' , 'J' ,...
%     'K' , 'L' , 'M' , 'N' , 'O' , 'P' , 'Q' , 'R' , 'S' , 'T' , 'U'...
%     , 'V' , 'W' , 'X' , 'Y' , 'Z'];
 capabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
 % Frequency vectors for finding the key once the length is determined.
 StdFreq = (1/100)*[8.167,1.492,2.782,4.253,12.702,2.228,2.015,6.094,6.966,...
     .0153,.0772,4.025,2.406,6.749,7.507,1.929,.095,5.987,6.327,...
     9.056,2.758,.978,2.36,.15,1.974,.074];
 Shift = zeros(26,26);
 % In this Shift matrix, each row corresponds to a possible key for a shift
 % cipher: the 3rd row represents the key where A -> D, B->E, etc.
 for i = 1:26
     Shift(i,:) = circshift(StdFreq,[0 i-1]);
 end
for i=1:26;
    message=strrep(message,capabet(i),alphabet(i));
end
X=1;
while X<=length(message),
    %This deletes all characters which are neither letters or spaces
    if (double(message(X))<97)&&not(double(message(X))==32);
        message(X)='';
        X=X+1;
    elseif (double(message(X))>122)&&not(double(message(X))==32);
        message(X)='';
        X=X+1;
    else X=X+1;
    end
end
keymax = input('Enter the maximum key length to check for: ');
    if keymax == ' '
        keymax = length(message);
    else
    end
    IOC = zeros(keymax,keymax);
% Now to check for various key lengths. The code will check key lengths
% from 1 (Caesar's Cipher) through user_input(keymax) and return the 
% diferent IOCs for each key length guess
for j = 1:keymax;
    submessage = (vec2mat(message,j));
    for kk=1:j;
        % Recording relative frequencies of each letter
            a = length(find(submessage(:,kk)=='a'));
            b = length(find(submessage(:,kk)=='b'));
            c = length(find(submessage(:,kk)=='c'));
            d = length(find(submessage(:,kk)=='d'));
            e = length(find(submessage(:,kk)=='e'));
            f = length(find(submessage(:,kk)=='f'));
            g = length(find(submessage(:,kk)=='g'));
            h = length(find(submessage(:,kk)=='h'));
            I = length(find(submessage(:,kk)=='i'));
            J = length(find(submessage(:,kk)=='j'));
            k = length(find(submessage(:,kk)=='k'));
            l = length(find(submessage(:,kk)=='l'));
            m = length(find(submessage(:,kk)=='n'));
            n = length(find(submessage(:,kk)=='m'));
            o = length(find(submessage(:,kk)=='o'));
            p = length(find(submessage(:,kk)=='p'));
            q = length(find(submessage(:,kk)=='q'));
            r = length(find(submessage(:,kk)=='r'));
            s = length(find(submessage(:,kk)=='s'));
            t = length(find(submessage(:,kk)=='t'));
            u = length(find(submessage(:,kk)=='u'));
            v = length(find(submessage(:,kk)=='v'));
            w = length(find(submessage(:,kk)=='w'));
            x = length(find(submessage(:,kk)=='x'));
            y = length(find(submessage(:,kk)=='y'));
            z = length(find(submessage(:,kk)=='z'));
        L=[a,b,c,d,e,f,g,h,I,J,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z];
        IOC(j,kk) = sum((L./length(submessage(:,kk))).^2)*26;
    end

end
    disp('Each row corresponds to the proposed key length.')
    disp('ie. the 4th row are the IOCs for guessing a key length of 4')
    disp(IOC)
    disp('IOC for English is around 1.76')
    keylength = input('Look at the table of IOCs. Enter the most probable key length: ');
    VigMessage = vec2mat(message,keylength); 
    CodeWord = (blanks(keylength));
    key = zeros(keylength, 26);
        for kk = 1:keylength
            % Recording relative frequencies of each letter in each column
            % after the message has been divided up into a matrix
            a = length(find(VigMessage(:,kk)=='a'));
            b = length(find(VigMessage(:,kk)=='b'));
            c = length(find(VigMessage(:,kk)=='c'));
            d = length(find(VigMessage(:,kk)=='d'));
            e = length(find(VigMessage(:,kk)=='e'));
            f = length(find(VigMessage(:,kk)=='f'));
            g = length(find(VigMessage(:,kk)=='g'));
            h = length(find(VigMessage(:,kk)=='h'));
            I = length(find(VigMessage(:,kk)=='i'));
            J = length(find(VigMessage(:,kk)=='j'));
            k = length(find(VigMessage(:,kk)=='k'));
            l = length(find(VigMessage(:,kk)=='l'));
            m = length(find(VigMessage(:,kk)=='n'));
            n = length(find(VigMessage(:,kk)=='m'));
            o = length(find(VigMessage(:,kk)=='o'));
            p = length(find(VigMessage(:,kk)=='p'));
            q = length(find(VigMessage(:,kk)=='q'));
            r = length(find(VigMessage(:,kk)=='r'));
            s = length(find(VigMessage(:,kk)=='s'));
            t = length(find(VigMessage(:,kk)=='t'));
            u = length(find(VigMessage(:,kk)=='u'));
            v = length(find(VigMessage(:,kk)=='v'));
            w = length(find(VigMessage(:,kk)=='w'));
            x = length(find(VigMessage(:,kk)=='x'));
            y = length(find(VigMessage(:,kk)=='y'));
            z = length(find(VigMessage(:,kk)=='z'));
            F=[a,b,c,d,e,f,g,h,I,J,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z];
                 DotProds = zeros(1,26);
                 for j = 1:26
                    DotProds(j) = dot(F,Shift(j,:));
                 end
                [~,idx] = max(DotProds(:));
                [X] = ind2sub(size(DotProds),idx);
                CodeWord(kk) = char('a'+X-1);
                key(kk,:) = Shift(j,:); 
        end
    TheCodeWord = ['The Code Word is: ', CodeWord];
    disp(TheCodeWord)
    Key = (CodeWord - 'a');
    disp(Key)
    for i = 1:keylength
        for j = 1:length(VigMessage(:,i))
            VigMessage(j,i) = char(mod(VigMessage(j,i)-97- Key(i),26)+97);
        end
    end
   attempt = reshape(VigMessage',[1,length(VigMessage)*keylength]);
   disp(attempt)
   toc
end
