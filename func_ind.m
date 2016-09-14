function indices = func_ind
% manually inserted "new-indices", as there was no pattern in the given data
% 1st column --> new indices      2nd column --> old indices

indices = [ 1    1
            2    2
            3    17
            4    38
            5    39
            6    40
            7    41
            8    18
            9    3
            10   4
            11   5
            12   6
            13   21
            14   32 
            15   26
            16   27
            17   33
            18   19
            19   7
            20   8
            21   9
            22   13
            23   20
            24   35
            25   28
            26   29
            27   34
            28   22
            29   14
            30   10
            31   11
            32   15
            33   23    
            34   36
            35   30
            36   31
            37   37
            38   24
            39   16
            40   12
            41   68
            42   69
            43   65
            44   44
            45   50
            46   45
            47   51
            48   66
            49   70
            50   71
            51   72
            52   73
            53   67
            54   52
            55   58
            56   59
            57   46
            58   62
            59   79
            60   74
            61   75
            62   80
            63   63
            64   47
            65   56
            66   57
            67   48
            68   60
            69   81
            70   76
            71   77
            72   82
            73   61
            74   53
            75   54
            76   49
            77   55
            78   64
            79   83
            80   78
            81   25
            82   43
            83   42
                    ];

%% Rearrange second time.. w/ more appropriate indexing

indices_new = indices_new2(indices(:,1));

indices = [indices_new(:,2) indices];
indices = sortrows(indices);

indices = [indices(:,1) indices(:,3)];

%%   Sanity Check

% s1 = sum(indices(:,1))
% s2 = sum(indices(:,2))
% p1 = prod(indices(:,1))
% p2 = prod(indices(:,2))
% 
% i1 = indices(:,1);
% i2 = indices(:,2);
% 
% u1 = unique(i1);
% u2 = unique(i2);
% 
% l1 = length(u1);
% l2 = length(u2);
% 
% diff = sort(i2)-i1;
