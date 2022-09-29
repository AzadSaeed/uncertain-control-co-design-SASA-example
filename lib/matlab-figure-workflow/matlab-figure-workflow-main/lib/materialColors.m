function C = materialColors

% create material color library (from https://materialui.co/colors)
C = createMaterialColors;

% visualize all the included colors
% greyflag = true; colorMatrix(C,greyflag)

end

% create material color library
function C = createMaterialColors

C.red = [
rgb(255, 235, 238)
rgb(255, 205, 210)
rgb(239, 154, 154)
rgb(229, 115, 115)
rgb(239, 83, 80)
rgb(244, 67, 54)
rgb(229, 57, 53)
rgb(211, 47, 47)
rgb(198, 40, 40)
rgb(183, 28, 28)
    ];

C.pink = [
rgb(252, 228, 236)
rgb(248, 187, 208)
rgb(244, 143, 177)
rgb(240, 98, 146)
rgb(236, 64, 122)
rgb(233, 30, 99)
rgb(216, 27, 96)
rgb(194, 24, 91)
rgb(173, 20, 87)
rgb(136, 14, 79)
    ];

C.purple = [
rgb(243, 229, 245)
rgb(225, 190, 231)
rgb(206, 147, 216)
rgb(186, 104, 200)
rgb(171, 71, 188)
rgb(156, 39, 176)
rgb(142, 36, 170)
rgb(123, 31, 162)
rgb(106, 27, 154)
rgb(74, 20, 140)
    ];

C.deeppurple = [
rgb(237, 231, 246)
rgb(209, 196, 233)
rgb(179, 157, 219)
rgb(149, 117, 205)
rgb(126, 87, 194)
rgb(103, 58, 183)
rgb(94, 53, 177)
rgb(81, 45, 168)
rgb(69, 39, 160)
rgb(49, 27, 146)
    ];

C.indigo = [
rgb(232, 234, 246)
rgb(197, 202, 233)
rgb(159, 168, 218)
rgb(121, 134, 203)
rgb(92, 107, 192)
rgb(63, 81, 181)
rgb(57, 73, 171)
rgb(48, 63, 159)
rgb(40, 53, 147)
rgb(26, 35, 126)
    ];

C.blue = [
rgb(227, 242, 253)
rgb(187, 222, 251)
rgb(144, 202, 249)
rgb(100, 181, 246)
rgb(66, 165, 245)
rgb(33, 150, 243)
rgb(30, 136, 229)
rgb(25, 118, 210)
rgb(21, 101, 192)
rgb(13, 71, 161)
    ];

C.lightblue = [
rgb(225, 245, 254)
rgb(179, 229, 252)
rgb(129, 212, 250)
rgb(79, 195, 247)
rgb(41, 182, 246)
rgb(3, 169, 244)
rgb(3, 155, 229)
rgb(2, 136, 209)
rgb(2, 119, 189)
rgb(1, 87, 155)
    ];

C.cyan = [
rgb(224, 247, 250)
rgb(178, 235, 242)
rgb(128, 222, 234)
rgb(77, 208, 225)
rgb(38, 198, 218)
rgb(0, 188, 212)
rgb(0, 172, 193)
rgb(0, 151, 167)
rgb(0, 131, 143)
rgb(0, 96, 100)
    ];

C.teal = [
rgb(224, 242, 241)
rgb(178, 223, 219)
rgb(128, 203, 196)
rgb(77, 182, 172)
rgb(38, 166, 154)
rgb(0, 150, 136)
rgb(0, 137, 123)
rgb(0, 121, 107)
rgb(0, 105, 92)
rgb(0, 77, 64)
    ];

C.green = [
rgb(232, 245, 233)
rgb(200, 230, 201)
rgb(165, 214, 167)
rgb(129, 199, 132)
rgb(102, 187, 106)
rgb(76, 175, 80)
rgb(67, 160, 71)
rgb(56, 142, 60)
rgb(46, 125, 50)
rgb(27, 94, 32)
    ];

C.lightgreen = [
rgb(241, 248, 233)
rgb(220, 237, 200)
rgb(197, 225, 165)
rgb(174, 213, 129)
rgb(156, 204, 101)
rgb(139, 195, 74)
rgb(124, 179, 66)
rgb(104, 159, 56)
rgb(85, 139, 47)
rgb(51, 105, 30)
    ];

C.lime = [
rgb(249, 251, 231)
rgb(240, 244, 195)
rgb(230, 238, 156)
rgb(220, 231, 117)
rgb(212, 225, 87)
rgb(205, 220, 57)
rgb(192, 202, 51)
rgb(175, 180, 43)
rgb(158, 157, 36)
rgb(130, 119, 23)
    ];

C.yellow = [
rgb(255, 253, 231)
rgb(255, 249, 196)
rgb(255, 245, 157)
rgb(255, 241, 118)
rgb(255, 238, 88)
rgb(255, 235, 59)
rgb(253, 216, 53)
rgb(251, 192, 45)
rgb(249, 168, 37)
rgb(245, 127, 23)
    ];

C.amber = [
rgb(255, 248, 225)
rgb(255, 236, 179)
rgb(255, 224, 130)
rgb(255, 213, 79)
rgb(255, 202, 40)
rgb(255, 193, 7)
rgb(255, 179, 0)
rgb(255, 160, 0)
rgb(255, 143, 0)
rgb(255, 111, 0)
    ];

C.orange = [
rgb(255, 243, 224)
rgb(255, 224, 178)
rgb(255, 204, 128)
rgb(255, 183, 77)
rgb(255, 167, 38)
rgb(255, 152, 0)
rgb(251, 140, 0)
rgb(245, 124, 0)
rgb(239, 108, 0)
rgb(230, 81, 0)
    ];

C.deeporange = [
rgb(251, 233, 231)
rgb(255, 204, 188)
rgb(255, 171, 145)
rgb(255, 138, 101)
rgb(255, 112, 67)
rgb(255, 87, 34)
rgb(244, 81, 30)
rgb(230, 74, 25)
rgb(216, 67, 21)
rgb(191, 54, 12)
    ];

C.brown = [
rgb(239, 235, 233)
rgb(215, 204, 200)
rgb(188, 170, 164)
rgb(161, 136, 127)
rgb(141, 110, 99)
rgb(121, 85, 72)
rgb(109, 76, 65)
rgb(93, 64, 55)
rgb(78, 52, 46)
rgb(62, 39, 35)
    ];

C.grey = [
rgb(250, 250, 250)
rgb(245, 245, 245)
rgb(238, 238, 238)
rgb(224, 224, 224)
rgb(189, 189, 189)
rgb(158, 158, 158)
rgb(117, 117, 117)
rgb(97, 97, 97)
rgb(66, 66, 66)
rgb(33, 33, 33)
    ];

C.bluegrey = [
rgb(236, 239, 241)
rgb(207, 216, 220)
rgb(176, 190, 197)
rgb(144, 164, 174)
rgb(120, 144, 156)
rgb(96, 125, 139)
rgb(84, 110, 122)
rgb(69, 90, 100)
rgb(55, 71, 79)
rgb(38, 50, 56)
    ];

end

% convert to 0-1 rgb vector format
function c = rgb(r,g,b)
c = [r g b]/255;
end

% convert rgb to greyscale
function c = rgb2grey(c)
grey = (0.299*c(1) + 0.587*c(2) + 0.114*c(3));
c = [grey grey grey];
end

% create matrix of all colors
function colorMatrix(C,greyflag)

% initialize figure
hf = figure; hf.Color= 'w'; hold on

% get field names
names = fieldnames(C);

% go through each field (color)s
for idx = 1:length(names)

    % extract
    c = C.(names{idx});

    % go through each color variant and plot
    for k = 1:length(c)
        if greyflag
            c_ = rgb2grey(c(k,:));
        else
            c_ = c(k,:);
        end
        rectangle('position',[idx-0.5 k-0.5 1 1],'FaceColor',c_);
    end

end

% customize axis
ha = gca;
ha.YTick = 1:k;
ha.XTick = 1:length(names);
ha.XTickLabel = names;
ha.YDir = 'reverse';
axis equal
xlim([0.5 idx+0.5])
ylim([0.5 length(c)+0.5])

end