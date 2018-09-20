function attempt=subdecode(message)
%This function uses frequency analysis to attack and decode messages
%encoded using a /substitution cipher/ which exchanges each letter of the
%alphabet for another. Unlike the famous "Caesar Shift" cipher, which
%"shifts" the letters of the alphabet along by a certain amount, giving
%just 26 possible ciphers, this has no specific order, giving 26!
%possibilities. The function is very unlikely to work first time but,
%depending of the length of the text, highly likely to get the letters "T"
%and "E" correct. From here the user may spot an error in the decoding,
%such as "TQE" which is obviously "THE", and the function offers them the
%chance to swap the offending letters as often as is required. After a few
%simple cases, the encrypted message is generally broken.
%It removes all punctuation and case from the original
%text.
 alphabet='abcdefghijklmnopqrstuvwxyz';
 capabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
if length(message)<50
    disp('Warning: message is short (less than 50 characters). Frequency analysis is unlikely to decode message. Attempt with a longer message.')
end
for i=1:26,
    message=strrep(message,capabet(i),alphabet(i));
end
X=1;
while X<=length(message),
    if (double(message(X))<97)&&not(double(message(X))==32);%This deletes all characters which are neither letters or spaces
        message(X)='';
        X=X+1;
    elseif (double(message(X))>122)&&not(double(message(X))==32);
        message(X)='';
        X=X+1;
    else X=X+1;
    end
end
%message has now been converted to all lower case wih no punctuation and
%may be attacked by frequency analysis
    a=0;b=0;c=0;d=0;e=0;f=0;g=0;h=0;I=0;J=0;k=0;l=0;m=0;n=0;o=0;p=0;q=0;r=0;s=0;t=0;u=0;v=0;w=0;x=0;y=0;z=0;
   
    for i=1:length(message),
        switch message(i) %recording relative frequencies of each letter
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
        F=length(mat2str(L))-53;
        analysis=cat(1,cat(2,'[a b c d e f g h i j k l m n o p q r s t u v w x y z]',zeros(1,F)),mat2str(L)) %displaying frequencies of each leter to user
        frequabet='ETAOINSHRDLUCMWFYGPBVKXJQZ'; %this is the alphabet organised by frequency of appearance in the english language
        
        %check for digrams before manipulating the code
        DigramLook = zeros(26,26);
        for k = 1:26
            for l = 1:26
                DigramFreq = length(strfind(message, strcat(char('A'+k-1), char('A'+l-1))));
                DigramLook(k,l) = DigramFreq;
            end
        end
        disp(DigramLook)
        
        for j=1:26
            C=max(L); %finds most frequent letter and assumes it must correspond to "e", then resets and repeats for "t", etc.
        for i=1:26,
            if L(i)==C,
                freq=i;
                message=strrep(message,char(freq+96),frequabet(j)); %replacing each letter in order of frequency with the most frequent letters in english
                L(i)=-1;
                break
            end
        end
        end
        disp('the most likely letters to be correct are "E" and "T".'),disp(' Look out for "T-E" which is likely the word "the", or an individual letter, most likely "A".'),disp(' If you find a potential mistake you may swap two letters. e.g TOE THH ->swap O and H->THE TOO') 
        message
        swap=input('would you like to swap two letters? (''''Y''''/''''N'''')');
        while swap=='Y', %offering user the chance to swap two letters
            lswap=input('Enter the two letters you would like to swap in the form [''''Letter1'''',''''Letter2'''']');
            message=strrep(message,lswap(1),'x');
            message=strrep(message,lswap(2),'y');
            message=strrep(message,'x',lswap(2));
            message=strrep(message,'y',lswap(1)) %displaying current state of decoded message
            swap=input('Would you like to swap another two letters? (''''Y''''/''''N'''')'); %offering a further chance to swap
        end
        attempt=message;
        end
        

            