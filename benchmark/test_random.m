function output = test_random(x)

    % 큰 행렬 생성

    A = rand(100, 100);

    

    % --- 메모리 체크 ---

    s = whos('A'); % A라는 변수의 정보를 가져옴

    fprintf('변수 A의 크기: %.2f MB\n', s.bytes / 1024^2);

    % -----------------

    

    output = x * A;

end
