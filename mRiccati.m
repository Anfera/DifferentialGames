function dXdt = mRiccati(t, P, A, S, Q, D)
P = reshape(P, size(Q)); %Convert from "n^2"-by-1 to "n"-by-"n"
dXdt = -D*P - P*A + P*S*P - Q; %Determine derivative
dXdt = dXdt(:); %Convert from "n"-by-"n" to "n^2"-by-1