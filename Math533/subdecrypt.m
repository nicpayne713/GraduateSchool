function attempt=subdecrypt(message)
%This function uses frequency analysis to attack and decode messages
%encoded using a /substitution cipher/ which exchanges each letter of the
%alphabet for another. Unlike the famous "Caesar Shift" cipher, which
%"shifts" the letters of the alphabet along by a certain amount, giving
%just 26 possible ciphers, this has no specific order, giving 26!
%possibilities. 
%It removes all punctuation and case from the original
%text.
global DigramLook DigramFreq 
 alphabet='abcdefghijklmnopqrstuvwxyz';
 Alphabet = ['A' , 'B' , 'C' , 'D' , 'E' , 'F' , 'G' , 'H', 'I' , 'J' ,...
     'K' , 'L' , 'M' , 'N' , 'O' , 'P' , 'Q' , 'R' , 'S' , 'T' , 'U'...
     , 'V' , 'W' , 'X' , 'Y' , 'Z'];
 capabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
if length(message)<50
  disp('Warning: message is short (less than 50 characters).Frequency analysis')
  disp('is unlikely to decode message. Attempt with a longer message.')
end
for i=1:26,
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
% Message has now been converted to all lower case wih no punctuation and
% May be attacked by frequency analysis
    a=0;b=0;c=0;d=0;e=0;f=0;g=0;h=0;I=0;J=0;k=0;l=0;m=0;n=0;o=0;
    p=0;q=0;r=0;s=0;t=0;u=0;v=0;w=0;x=0;y=0;z=0;
   
    for i=1:length(message),
        % Recording relative frequencies of each letter
        switch message(i) 
            case 'a'
                a=a+1;
            case 'b'
                b=b+1;
            case 'c'
                c=c+1;
            case 'd'
                d=d+1;
            case 'e'
                e=e+1;
            case 'f'
                f=f+1;
            case 'g'
                g=g+1;
            case 'h'
                h=h+1;
            case 'i'
                I=I+1;
            case 'j'
                J=J+1;
            case 'k'
                k=k+1;
            case 'l'
                l=l+1;
            case 'm'
                m=m+1;
            case 'n'
                n=n+1;
            case 'o'
                o=o+1;
            case 'p'
                p=p+1;
            case 'q'
                q=q+1;
            case 'r'
                r=r+1;
            case 's'
                s=s+1;
            case 't'
                t=t+1;
            case 'u'
                u=u+1;
            case 'v'
                v=v+1;
            case 'w'
                w=w+1;
            case 'x'
                x=x+1;
            case 'y'
                y=y+1;
            case 'z'
                z=z+1;
        end
    end
        L=[a,b,c,d,e,f,g,h,I,J,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z];
        Character = Alphabet';
        Frequency = L';
        Frequencies = table(Character, Frequency);
        disp(Frequencies)
        
        % Calculate the IOC
        IOC = sum((L./length(message)).^2)*26; 
        disp('IOC = ') ,disp(IOC) 

        % Check for digrams before manipulating the code
        DigramLook = zeros(26,26);
        for k = 1:26
            for l = 1:26
                DigramFreq = length(strfind(message,...
                    strcat(char('a'+k-1), char('a'+l-1))));
                DigramLook(k,l) = DigramFreq;
            end
        end
        DigramMin = input('Input the minimum frequency a digram appears to display: ');
        C = find(DigramLook > DigramMin);
        DigramFreqs = zeros(1,length(C));
        Digrams = num2str(zeros(length(C),2));
        for i = 1:length(C)
            [q,r] = quorem(sym(C(i)),sym(26));
            Digrams(i,1) = char('A'+double(r-1));
            Digrams(i,2) = char('A'+double(q));
            DigramFreqs(i) = (DigramLook(double(r),double(q+1)));
        end
        
        DigramFreqs = DigramFreqs';
        CommonDigramsInCipherText = table(Digrams,DigramFreqs);
        disp(CommonDigramsInCipherText)      

attempt = message;

swap=input('would you like to swap letters? (''''Y''''/''''N'''')');
    while swap=='Y', %offering user the chance to swap two letters  
    disp('Enter the letter from the message in LOWERCASE then the letter you think it gets')
    lswap=input('sent to in the form [''Letter1'',''Letter2'';''Letter3'',''Letter4''; etc.]'...
            );
        for i = 1:length(message);
            for j = 1:length(lswap')
                if message(i) == lswap(j,1);
                      attempt(i) = lswap(j,2);
                end
            end
        end
    disp(attempt)
    swap=input('would you like to swap letters? (''''Y''''/''''N'''')');
    end
disp(attempt)    
end

        
        
        
        
        
        
        