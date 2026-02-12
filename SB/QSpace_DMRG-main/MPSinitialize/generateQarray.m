function arr = generateQarray(Nsite, S)
    % Check if S is reachable
    if S > Nsite || S < 0
        error('S must be between 0 and Nsite');
    end

    if mod(Nsite+S,2) ~=0 % Parity Check
        error('S and Nsite must have same parity');
    end
    
    % Initialize array
    arr = zeros(1, Nsite + 1);
    arr(1) = 0;
    
    % Calculate the number of increases and decreases
    increases = (Nsite + S)/2;
    decreases = (Nsite-S)/2;
    
    % Shuffle the sequence
    sequence = generateValidPermutation(increases, decreases,Nsite,false);
    
    % Apply the sequence while ensuring constraints
    for i = 2:Nsite+1
        arr(i) = arr(i-1)+sequence(i-1);
    end
end

function seq = generateValidPermutation(inc, dec,N,finIs1)
    if inc-dec < 0 || inc-dec > N
        seq = NaN;
        return;
    end
    if N == 1
        seq = 1;
        return;
    end
    if finIs1
        fin = 1;
    else
        fin = randomChoice();
    end
    seq = [generateValidPermutation(inc-(fin == 1),dec-(fin == -1),N-1,(fin==-1)) fin];
    if isnan(seq(1))
        fin  = -fin;
        seq = [generateValidPermutation(inc-(fin == 1),dec-(fin == -1),N-1,(fin==-1)) fin];
    end
end

function choice = randomChoice()
    if rand() < 0.5
        choice = 1;
    else
        choice = -1;
    end
end