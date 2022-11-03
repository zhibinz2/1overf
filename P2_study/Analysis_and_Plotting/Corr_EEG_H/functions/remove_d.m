function [Y,H]=remove_d(X,d)

    ftransf=reshape(fft(X),[],1);
    [H1,~]=freqz([1,-1],1,length(X),'whole');
    H=H1.^d;
    u=H.*ftransf;
    if d<0, u=[u(2); u(2:length(u))]; end %Brockwell and Davis 1981 pag 530
    Y=real(ifft(u));
    
end