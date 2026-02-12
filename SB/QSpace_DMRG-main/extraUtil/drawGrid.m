function h = drawGrid(Xarray, Yarray, data, varargin)
    % drawGrid - 데이터를 그리드 형태로 표시하는 함수
    %
    % 사용법:
    %   h = drawGrid(Xarray, Yarray, data)
    %   h = drawGrid(Xarray, Yarray, data, 'Name', Value, ...)
    %
    % 입력:
    %   Xarray - X 축 값 배열
    %   Yarray - Y 축 값 배열
    %   data   - 그리드 데이터
    %
    % 선택적 입력 (이름-값 쌍):
    %   'Colormap'  - 컬러맵 (기본값: parula)
    %   'CLim'      - 컬러 스케일 범위 [min max] (기본값: 자동)
    %   'NanColor'  - NaN 값의 색상 [R G B] (기본값: [1 1 1] 흰색)
    %   'NanWidth'  - NaN 색상의 두께 (기본값: 0.02, 즉 컬러맵의 2%)
    
    % 입력 파싱
    p = inputParser;
    addRequired(p, 'Xarray');
    addRequired(p, 'Yarray');
    addRequired(p, 'data');
    addParameter(p, 'Colormap', parula(256));
    addParameter(p, 'CLim', []);
    addParameter(p, 'NanColor', [1 1 1]);
    addParameter(p, 'NanWidth', 0.02); % NaN 색상 두께 (컬러맵의 비율)
    parse(p, Xarray, Yarray, data, varargin{:});
    
    % 파싱된 결과 가져오기
    cmap = p.Results.Colormap;
    clim = p.Results.CLim;
    nanColor = p.Results.NanColor;
    nanWidth = p.Results.NanWidth;
    
    % NaN이 아닌 데이터의 최대값과 최소값 찾기
    valid_data = data(~isnan(data));
    if isempty(clim)
        clim = [min(valid_data(:)) max(valid_data(:))];
    end
    
    % 데이터 복사
    data_masked = data;
    
    % 컬러맵 크기 가져오기
    n_colors = size(cmap, 1);
    
    % NaN 색상 두께를 컬러맵 인덱스 개수로 변환
    nan_size = max(1, round(n_colors * nanWidth));
    
    % 새로운 컬러맵 생성 (NaN 색상을 맨 위에 추가)
    new_cmap = [cmap; repmat(nanColor, nan_size, 1)];
    
    % NaN 값을 컬러맵의 맨 위 값으로 대체
    nan_mask = isnan(data);
    data_masked(nan_mask) = clim(2) + 1; % 최대값보다 약간 더 큰 값 사용
    
    % 이미지 표시
    h = imagesc(Xarray, Yarray, data_masked);
    colormap(gca, new_cmap);
    
    % 컬러바 범위 확장 (NaN 값을 포함하도록)
    % NaN 영역이 전체 컬러맵에서 차지하는 비율 계산
    nan_ratio = nan_size / (n_colors + nan_size);
    
    % 확장된 범위 계산
    extended_range = (clim(2) - clim(1)) / (1 - nan_ratio);
    new_max = clim(1) + extended_range;
    
    % 컬러바 범위 설정
    caxis([clim(1), new_max]);
    
    % 컬러바 표시 및 눈금 조정
    c = colorbar;
    
    % 컬러바 눈금 설정 (NaN 영역 포함)
    n_ticks = 5;
    tick_values = linspace(clim(1), clim(2), n_ticks);
    tick_positions = linspace(clim(1), clim(1) + (1-nan_ratio)*extended_range, n_ticks);
    
    % NaN 눈금 위치 (증가하는 순서로 추가해야 함)
    nan_tick_position = clim(1) + extended_range * 0.98; % NaN 영역 중간쯤
    
    % 모든 tick 위치를 정렬된 배열로 만들기
    all_positions = sort([tick_positions, nan_tick_position]);
    
    % 레이블 생성 (정렬된 tick 위치에 맞춰서)
    all_labels = cell(size(all_positions));
    for i = 1:length(all_positions)
        if all_positions(i) == nan_tick_position
            all_labels{i} = 'NaN';
        else
            % 가장 가까운 원래 tick 값 찾기
            [~, idx] = min(abs(tick_positions - all_positions(i)));
            all_labels{i} = sprintf('%.2f', tick_values(idx));
        end
    end
    
    c.Ticks = all_positions;
    c.TickLabels = all_labels;
end