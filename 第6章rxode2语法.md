# 第六章rxode2语法

本章节简要描述了用于定义模型的语法，`rxode2`将把这些模型转换为r可调用的编译代码。它还描述了`R`和`rxode2`模型指定之间变量的通信。

## 6.1示例

```R
   # An rxode2 model specification (this line is a comment).

   if(comed==0){   # concomitant medication (con-med)?
      F = 1.0;     # full bioavailability w.o. con-med
   }
   else {
      F = 0.80;    # 20% reduced bioavailability
   }

   C2 = centr/V2;  # concentration in the central compartment
   C3 = peri/V3;   # concentration in the peripheral compartment

   # ODE describing the PK and PD

   d/dt(depot) = -KA*depot;
   d/dt(centr) = F*KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                      Q*C2 - Q*C3;
   d/dt(eff)   = Kin - Kout*(1-C2/(EC50+C2))*eff;
```

## 6.2语法

一个`rxode2`模型描述由一个或多个语句组成，不同语句间可在句尾可选地添加分号`;`进行分隔，和可选添加注释(注释由#开始和在行尾结束)。

语句块是由大括号`{...}`分隔的一组语句。

语句可以是赋值、条件`if`/`else if`/`else`、`while`循环(可以通过`break`退出)、特殊语句或打印语句(用于调试/测试)。

赋值语句可以是：

- **简单**赋值，其中赋值符左手边是标识符(即变量)
- 对特殊的**时间导数**赋值，其中赋值符左手边指定相应状态变量(房室)中的数量相对于时间的变化，例如，`d/dt(depot)`：
- 对特殊的**初始条件**赋值，其中赋值符左手边指定被指定的初始条件的房室， 例如`depot(0) = 0`
- 对特殊的模型事件的修改，包括**生物利用度** (`f(depot)=1`)，**滞后时间**(`alag(depot)=0`)，**使用输注速率参数化** (`rate(depot)=2`)和**使用输注持续时长参数化**(`dur(depot)=2`)。在[rxode2事件部分](https://nlmixr2.github.io/rxode2-manual/events.html#events)中可以找到这些模型特征的一个示例，以及用于描述输注参数化方式的事件规范和rxode2数据规范。
- 特殊**更改点语法或模型时间**。这些模型时间是由`mtime(var)=time`指定
- 特殊的**雅可比导数**复制，其中赋值符左手边指定房室ode相相对于某个变量的变化。例如，如果`d/dt(y) = dy`，则为此使用Jacobian 房室可以指定为`df(y)/dy(dy) = 1`。对于非常刚性的ODE系统，获得解或指定雅可比矩阵可能会有一些好处。然而，对于我们使用LSODA尝试的几个刚性系统，这实际上略微减慢了求解速度。

注意赋值可以通过`=`、`<-`或`~`来完成。

使用`~`运算符赋值时，不会输出**简单赋值**和**时间导数**赋值。

特别声明可以是：

- **房室声明语句**，它可以改变默认的给药房室和假设的房室数，以及在最后添加额外的房室名(对多终点nlmixr模型有用);这些是由`cmt(compartmentName)`指定
- **参数声明语句**, 它可以确保输入参数按一定的顺序排列，而不是按解析的顺序排列参数。这对于在使用2个不同的ODE模型时保持参数顺序相同非常有用。它们由 `param(par1, par2,...)`指定

一个示例模型如下所示:

```R
   # simple assignment
   C2 = centr/V2;

   # time-derivative assignment
   d/dt(centr) = F*KA*depot - CL*C2 - Q*C2 + Q*C3; 
```

赋值和`if`语句中的表达式可以是数值或逻辑值。

数字表达式可以包括以下数字运算符`+, -, *, /, ^`以及C或R数学库中定义的数学函数(例如`fabs`、`exp`、`log`、`sin`、`abs`)。

您还可以访问[R数学库](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Numerical-analysis-subroutines)中的R函数， 就像log gamma函数的`lgammafn`一样。

`rxode2`语法区分大小写，即`ABC`与`abc`、`Abc`、`ABc`等是不同。

### 6.2.1 标识符

与R一样，标识符(变量名)可能由一个或多个字母数字， 下划线`_`或句点`.`字符组成，但第一个字符不能是数字或下划线`_`。

模型描述中的标识符可以参考：

- 动态系统中的状态变量(例如，药代动力学模型中的房室)。
- 隐含输入变量`t`(时间)、`tlast`(最后时间点)和 `podo`(口服剂量，在未记录的吸收过转移室型中)。
- 特殊常数如`pi`或 [R的预定义常量](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Mathematical-constants)。
- 模型参数(例如`ka`吸收率、`CL`清除率等)
- 其他，作为模型描述的一部分由赋值符创建的;这些被称为LHS(左手边)变量。

目前`rxode2`建模语言只识别系统状态变量和“参数”，因此，需要从R传递到ODE模型的任何值(例如，`age`)都应该在`rxSolve()`的`params`参数中传递，或者在提供的事件数据集中传递。

在`rxode2`事件表中有某些变量名。 为避免混淆，以下事件表相关项目不能 被赋值或用作状态，但可以在rxode2代码中访问：

- `cmt`
- `dvid`
- `addl`
- `ss`
- `rate`
- `id`

但是，以下变量不能在模型描述中使用：

- `evid`
- `ii`

有时rxode2会生成反馈给rxode2的变量。 类似地，nlmixr生成一些变量，用于nlmixr估计和模拟。这些变量以 `rx`或`nlmixr`前缀。为避免任何问题，建议不要为变量使用`rx`或`nlmixr` 前缀。

## 6.3逻辑运算符

逻辑运算符支持标准的R运算符`==`、`!=` `>=` `<=` `>`和`<`。像R一样，它们可以用在`if()`或`while()` 语句，`ifelse()`表达式。此外，它们可以在标准的赋值中。例如，以下内容有效：

```R
cov1 = covm*(sexf == "female") + covm*(sexf != "female")
```

请注意，您还可以在比较中使用字符表达式。 这种便利是有代价的，因为字符比较比数值表达式慢。与R不同，`as.numeric`或 `as.integer`对于这些逻辑语句不仅不需要，而且如果您尝试使用这些函数，将导致语法错误。

## 6.4cmt()改变状态的房室编号

使用`cmt()`语句可以更改房室的顺序。要了解`cmt()`可以做什么，您需要了解`rxode2`如何编号房室。

下面是rxode2如何编号房室的示例

### 6.4.1rxode2如何编号房室

rxode2在解析时自动分配房室编号。例如，对于Mavoglurant PBPK模型，可以使用以下模型：

```R
library(rxode2)
pbpk <- rxode2({
    KbBR = exp(lKbBR)
    KbMU = exp(lKbMU)
    KbAD = exp(lKbAD)
    CLint= exp(lCLint + eta.LClint)
    KbBO = exp(lKbBO)
    KbRB = exp(lKbRB)

    ## Regional blood flows
    # Cardiac output (L/h) from White et al (1968)
    CO  = (187.00*WT^0.81)*60/1000; 
    QHT = 4.0 *CO/100;
    QBR = 12.0*CO/100;
    QMU = 17.0*CO/100;
    QAD = 5.0 *CO/100;
    QSK = 5.0 *CO/100;
    QSP = 3.0 *CO/100;
    QPA = 1.0 *CO/100;
    QLI = 25.5*CO/100;
    QST = 1.0 *CO/100;
    QGU = 14.0*CO/100;
    # Hepatic artery blood flow
    QHA = QLI - (QSP + QPA + QST + QGU); 
    QBO = 5.0 *CO/100;
    QKI = 19.0*CO/100;
    QRB = CO - (QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI);
    QLU = QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI + QRB;

    ## Organs' volumes = organs' weights / organs' density
    VLU = (0.76 *WT/100)/1.051;
    VHT = (0.47 *WT/100)/1.030;
    VBR = (2.00 *WT/100)/1.036;
    VMU = (40.00*WT/100)/1.041;
    VAD = (21.42*WT/100)/0.916;
    VSK = (3.71 *WT/100)/1.116;
    VSP = (0.26 *WT/100)/1.054;
    VPA = (0.14 *WT/100)/1.045;
    VLI = (2.57 *WT/100)/1.040;
    VST = (0.21 *WT/100)/1.050;
    VGU = (1.44 *WT/100)/1.043;
    VBO = (14.29*WT/100)/1.990;
    VKI = (0.44 *WT/100)/1.050;
    VAB = (2.81 *WT/100)/1.040;
    VVB = (5.62 *WT/100)/1.040;
    VRB = (3.86 *WT/100)/1.040;

    ## Fixed parameters
    BP = 0.61;      # Blood:plasma partition coefficient
    fup = 0.028;    # Fraction unbound in plasma
    fub = fup/BP;   # Fraction unbound in blood

    KbLU = exp(0.8334);
    KbHT = exp(1.1205);
    KbSK = exp(-.5238);
    KbSP = exp(0.3224);
    KbPA = exp(0.3224);
    KbLI = exp(1.7604);
    KbST = exp(0.3224);
    KbGU = exp(1.2026);
    KbKI = exp(1.3171);

    ##-----------------------------------------
    S15 = VVB*BP/1000;
    C15 = Venous_Blood/S15

    ##-----------------------------------------
    d/dt(Lungs) = QLU*(Venous_Blood/VVB - Lungs/KbLU/VLU);
    d/dt(Heart) = QHT*(Arterial_Blood/VAB - Heart/KbHT/VHT);
    d/dt(Brain) = QBR*(Arterial_Blood/VAB - Brain/KbBR/VBR);
    d/dt(Muscles) = QMU*(Arterial_Blood/VAB - Muscles/KbMU/VMU);
    d/dt(Adipose) = QAD*(Arterial_Blood/VAB - Adipose/KbAD/VAD);
    d/dt(Skin) = QSK*(Arterial_Blood/VAB - Skin/KbSK/VSK);
    d/dt(Spleen) = QSP*(Arterial_Blood/VAB - Spleen/KbSP/VSP);
    d/dt(Pancreas) = QPA*(Arterial_Blood/VAB - Pancreas/KbPA/VPA);
    d/dt(Liver) = QHA*Arterial_Blood/VAB + QSP*Spleen/KbSP/VSP +
      QPA*Pancreas/KbPA/VPA + QST*Stomach/KbST/VST +
      QGU*Gut/KbGU/VGU - CLint*fub*Liver/KbLI/VLI - QLI*Liver/KbLI/VLI;
    d/dt(Stomach) = QST*(Arterial_Blood/VAB - Stomach/KbST/VST);
    d/dt(Gut) = QGU*(Arterial_Blood/VAB - Gut/KbGU/VGU);
    d/dt(Bones) = QBO*(Arterial_Blood/VAB - Bones/KbBO/VBO);
    d/dt(Kidneys) = QKI*(Arterial_Blood/VAB - Kidneys/KbKI/VKI);
    d/dt(Arterial_Blood) = QLU*(Lungs/KbLU/VLU - Arterial_Blood/VAB);
    d/dt(Venous_Blood) = QHT*Heart/KbHT/VHT + QBR*Brain/KbBR/VBR +
      QMU*Muscles/KbMU/VMU + QAD*Adipose/KbAD/VAD + QSK*Skin/KbSK/VSK +
      QLI*Liver/KbLI/VLI + QBO*Bones/KbBO/VBO + QKI*Kidneys/KbKI/VKI +
      QRB*Rest_of_Body/KbRB/VRB - QLU*Venous_Blood/VVB;
    d/dt(Rest_of_Body) = QRB*(Arterial_Blood/VAB - Rest_of_Body/KbRB/VRB);
})
```

如果您查看摘要，您可以看到rxode2赋值房室号的位置

```R
summary(pbpk)
```

```R
#> rxode2 2.0.11 model named rx_291007fc063b6ec76a6a1e59198481c3 model (✔ ready). 
#> DLL: /home/matt/.cache/R/rxode2/rx_291007fc063b6ec76a6a1e59198481c3__.rxd/rx_291007fc063b6ec76a6a1e59198481c3_.so
#> NULL
#> 
#> Calculated Variables:
#>  [1] "KbBR"  "KbMU"  "KbAD"  "CLint" "KbBO"  "KbRB"  "CO"    "QHT"   "QBR"  
#> [10] "QMU"   "QAD"   "QSK"   "QSP"   "QPA"   "QLI"   "QST"   "QGU"   "QHA"  
#> [19] "QBO"   "QKI"   "QRB"   "QLU"   "VLU"   "VHT"   "VBR"   "VMU"   "VAD"  
#> [28] "VSK"   "VSP"   "VPA"   "VLI"   "VST"   "VGU"   "VBO"   "VKI"   "VAB"  
#> [37] "VVB"   "VRB"   "fub"   "KbLU"  "KbHT"  "KbSK"  "KbSP"  "KbPA"  "KbLI" 
#> [46] "KbST"  "KbGU"  "KbKI"  "S15"   "C15"  
#> ── rxode2 Model Syntax ──
#> rxode2({
#>     KbBR = exp(lKbBR)
#>     KbMU = exp(lKbMU)
#>     KbAD = exp(lKbAD)
#>     CLint = exp(lCLint + eta.LClint)
#>     KbBO = exp(lKbBO)
#>     KbRB = exp(lKbRB)
#>     CO = (187 * WT^0.81) * 60/1000
#>     QHT = 4 * CO/100
#>     QBR = 12 * CO/100
#>     QMU = 17 * CO/100
#>     QAD = 5 * CO/100
#>     QSK = 5 * CO/100
#>     QSP = 3 * CO/100
#>     QPA = 1 * CO/100
#>     QLI = 25.5 * CO/100
#>     QST = 1 * CO/100
#>     QGU = 14 * CO/100
#>     QHA = QLI - (QSP + QPA + QST + QGU)
#>     QBO = 5 * CO/100
#>     QKI = 19 * CO/100
#>     QRB = CO - (QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI)
#>     QLU = QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI + QRB
#>     VLU = (0.76 * WT/100)/1.051
#>     VHT = (0.47 * WT/100)/1.03
#>     VBR = (2 * WT/100)/1.036
#>     VMU = (40 * WT/100)/1.041
#>     VAD = (21.42 * WT/100)/0.916
#>     VSK = (3.71 * WT/100)/1.116
#>     VSP = (0.26 * WT/100)/1.054
#>     VPA = (0.14 * WT/100)/1.045
#>     VLI = (2.57 * WT/100)/1.04
#>     VST = (0.21 * WT/100)/1.05
#>     VGU = (1.44 * WT/100)/1.043
#>     VBO = (14.29 * WT/100)/1.99
#>     VKI = (0.44 * WT/100)/1.05
#>     VAB = (2.81 * WT/100)/1.04
#>     VVB = (5.62 * WT/100)/1.04
#>     VRB = (3.86 * WT/100)/1.04
#>     BP = 0.61
#>     fup = 0.028
#>     fub = fup/BP
#>     KbLU = exp(0.8334)
#>     KbHT = exp(1.1205)
#>     KbSK = exp(-0.5238)
#>     KbSP = exp(0.3224)
#>     KbPA = exp(0.3224)
#>     KbLI = exp(1.7604)
#>     KbST = exp(0.3224)
#>     KbGU = exp(1.2026)
#>     KbKI = exp(1.3171)
#>     S15 = VVB * BP/1000
#>     C15 = Venous_Blood/S15
#>     d/dt(Lungs) = QLU * (Venous_Blood/VVB - Lungs/KbLU/VLU)
#>     d/dt(Heart) = QHT * (Arterial_Blood/VAB - Heart/KbHT/VHT)
#>     d/dt(Brain) = QBR * (Arterial_Blood/VAB - Brain/KbBR/VBR)
#>     d/dt(Muscles) = QMU * (Arterial_Blood/VAB - Muscles/KbMU/VMU)
#>     d/dt(Adipose) = QAD * (Arterial_Blood/VAB - Adipose/KbAD/VAD)
#>     d/dt(Skin) = QSK * (Arterial_Blood/VAB - Skin/KbSK/VSK)
#>     d/dt(Spleen) = QSP * (Arterial_Blood/VAB - Spleen/KbSP/VSP)
#>     d/dt(Pancreas) = QPA * (Arterial_Blood/VAB - Pancreas/KbPA/VPA)
#>     d/dt(Liver) = QHA * Arterial_Blood/VAB + QSP * Spleen/KbSP/VSP + 
#>         QPA * Pancreas/KbPA/VPA + QST * Stomach/KbST/VST + QGU * 
#>         Gut/KbGU/VGU - CLint * fub * Liver/KbLI/VLI - QLI * Liver/KbLI/VLI
#>     d/dt(Stomach) = QST * (Arterial_Blood/VAB - Stomach/KbST/VST)
#>     d/dt(Gut) = QGU * (Arterial_Blood/VAB - Gut/KbGU/VGU)
#>     d/dt(Bones) = QBO * (Arterial_Blood/VAB - Bones/KbBO/VBO)
#>     d/dt(Kidneys) = QKI * (Arterial_Blood/VAB - Kidneys/KbKI/VKI)
#>     d/dt(Arterial_Blood) = QLU * (Lungs/KbLU/VLU - Arterial_Blood/VAB)
#>     d/dt(Venous_Blood) = QHT * Heart/KbHT/VHT + QBR * Brain/KbBR/VBR + 
#>         QMU * Muscles/KbMU/VMU + QAD * Adipose/KbAD/VAD + QSK * 
#>         Skin/KbSK/VSK + QLI * Liver/KbLI/VLI + QBO * Bones/KbBO/VBO + 
#>         QKI * Kidneys/KbKI/VKI + QRB * Rest_of_Body/KbRB/VRB - 
#>         QLU * Venous_Blood/VVB
#>     d/dt(Rest_of_Body) = QRB * (Arterial_Blood/VAB - Rest_of_Body/KbRB/VRB)
#> })
```

在此案例中，`Venous_Blood`被分配到房室`15`。 弄清楚这一点可能会很不方便，并且还会导致模拟或估计数据集中的房室重新编号。虽然通过名称指定房室很容易而且可能更清晰，但其他工具只支持房室编号。因此，有一种方法可以很容易地为房室编号，从而减少多个工具之间的数据修改。

### 6.4.2 使用`cmt()`预声明更改房室

要按您希望的顺序将房室添加到rxode2模型中，您只需要使用cmt预先声明这些房室。例如，将Venous_Blood和Skin分别指定为第一个和第二个房室是很简单的:

```R
pbpk2 <- rxode2({
  ## Now this is the first compartment, ie cmt=1
  cmt(Venous_Blood)
  ## Skin may be a compartment you wish to dose to as well,
  ##  so it is now cmt=2
  cmt(Skin) 
  KbBR = exp(lKbBR)
  KbMU = exp(lKbMU)
  KbAD = exp(lKbAD)
  CLint= exp(lCLint + eta.LClint)
  KbBO = exp(lKbBO)
  KbRB = exp(lKbRB)

  ## Regional blood flows
  # Cardiac output (L/h) from White et al (1968)m
  CO  = (187.00*WT^0.81)*60/1000; 
  QHT = 4.0 *CO/100;
  QBR = 12.0*CO/100;
  QMU = 17.0*CO/100;
  QAD = 5.0 *CO/100;
  QSK = 5.0 *CO/100;
  QSP = 3.0 *CO/100;
  QPA = 1.0 *CO/100;
  QLI = 25.5*CO/100;
  QST = 1.0 *CO/100;
  QGU = 14.0*CO/100;
  QHA = QLI - (QSP + QPA + QST + QGU); # Hepatic artery blood flow
  QBO = 5.0 *CO/100;
  QKI = 19.0*CO/100;
  QRB = CO - (QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI);
  QLU = QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI + QRB;

  ## Organs' volumes = organs' weights / organs' density
  VLU = (0.76 *WT/100)/1.051;
  VHT = (0.47 *WT/100)/1.030;
  VBR = (2.00 *WT/100)/1.036;
  VMU = (40.00*WT/100)/1.041;
  VAD = (21.42*WT/100)/0.916;
  VSK = (3.71 *WT/100)/1.116;
  VSP = (0.26 *WT/100)/1.054;
  VPA = (0.14 *WT/100)/1.045;
  VLI = (2.57 *WT/100)/1.040;
  VST = (0.21 *WT/100)/1.050;
  VGU = (1.44 *WT/100)/1.043;
  VBO = (14.29*WT/100)/1.990;
  VKI = (0.44 *WT/100)/1.050;
  VAB = (2.81 *WT/100)/1.040;
  VVB = (5.62 *WT/100)/1.040;
  VRB = (3.86 *WT/100)/1.040;

  ## Fixed parameters
  BP = 0.61;      # Blood:plasma partition coefficient
  fup = 0.028;    # Fraction unbound in plasma
  fub = fup/BP;   # Fraction unbound in blood

  KbLU = exp(0.8334);
  KbHT = exp(1.1205);
  KbSK = exp(-.5238);
  KbSP = exp(0.3224);
  KbPA = exp(0.3224);
  KbLI = exp(1.7604);
  KbST = exp(0.3224);
  KbGU = exp(1.2026);
  KbKI = exp(1.3171);


  ##-----------------------------------------
  S15 = VVB*BP/1000;
  C15 = Venous_Blood/S15

  ##-----------------------------------------
  d/dt(Lungs) = QLU*(Venous_Blood/VVB - Lungs/KbLU/VLU);
  d/dt(Heart) = QHT*(Arterial_Blood/VAB - Heart/KbHT/VHT);
  d/dt(Brain) = QBR*(Arterial_Blood/VAB - Brain/KbBR/VBR);
  d/dt(Muscles) = QMU*(Arterial_Blood/VAB - Muscles/KbMU/VMU);
  d/dt(Adipose) = QAD*(Arterial_Blood/VAB - Adipose/KbAD/VAD);
  d/dt(Skin) = QSK*(Arterial_Blood/VAB - Skin/KbSK/VSK);
  d/dt(Spleen) = QSP*(Arterial_Blood/VAB - Spleen/KbSP/VSP);
  d/dt(Pancreas) = QPA*(Arterial_Blood/VAB - Pancreas/KbPA/VPA);
  d/dt(Liver) = QHA*Arterial_Blood/VAB + QSP*Spleen/KbSP/VSP +
    QPA*Pancreas/KbPA/VPA + QST*Stomach/KbST/VST + QGU*Gut/KbGU/VGU -
    CLint*fub*Liver/KbLI/VLI - QLI*Liver/KbLI/VLI;
  d/dt(Stomach) = QST*(Arterial_Blood/VAB - Stomach/KbST/VST);
  d/dt(Gut) = QGU*(Arterial_Blood/VAB - Gut/KbGU/VGU);
  d/dt(Bones) = QBO*(Arterial_Blood/VAB - Bones/KbBO/VBO);
  d/dt(Kidneys) = QKI*(Arterial_Blood/VAB - Kidneys/KbKI/VKI);
  d/dt(Arterial_Blood) = QLU*(Lungs/KbLU/VLU - Arterial_Blood/VAB);
  d/dt(Venous_Blood) = QHT*Heart/KbHT/VHT + QBR*Brain/KbBR/VBR +
    QMU*Muscles/KbMU/VMU + QAD*Adipose/KbAD/VAD + QSK*Skin/KbSK/VSK +
    QLI*Liver/KbLI/VLI + QBO*Bones/KbBO/VBO + QKI*Kidneys/KbKI/VKI +
    QRB*Rest_of_Body/KbRB/VRB - QLU*Venous_Blood/VVB;
  d/dt(Rest_of_Body) = QRB*(Arterial_Blood/VAB - Rest_of_Body/KbRB/VRB);
})
```

您可以在简单的打印输出中看到此更改

```R
pbpk2
```

```R
#> rxode2 2.0.11 model named rx_8538903f734422ef88399de66a046870 model (✔ ready). 
#> x$state: Venous_Blood, Skin, Lungs, Heart, Brain, Muscles, Adipose, Spleen, Pancreas, Liver, Stomach, Gut, Bones, Kidneys, Arterial_Blood, Rest_of_Body
#> x$params: lKbBR, lKbMU, lKbAD, lCLint, eta.LClint, lKbBO, lKbRB, WT, BP, fup
#> x$lhs: KbBR, KbMU, KbAD, CLint, KbBO, KbRB, CO, QHT, QBR, QMU, QAD, QSK, QSP, QPA, QLI, QST, QGU, QHA, QBO, QKI, QRB, QLU, VLU, VHT, VBR, VMU, VAD, VSK, VSP, VPA, VLI, VST, VGU, VBO, VKI, VAB, VVB, VRB, fub, KbLU, KbHT, KbSK, KbSP, KbPA, KbLI, KbST, KbGU, KbKI, S15, C15
```

前两个房室是`Venous_Blood`后面`Skin`。

### 6.4.3用`cmt()`将房室附加到模型

您还可以将“房室”附加到模型中。由于ODE求解内部结构，在定义所有微分方程之前，您不能将假的房室添加到模型中。

例如，这是合法的：

```R
ode.1c.ka <- rxode2({
    C2 = center/V;
    d / dt(depot) = -KA * depot
    d/dt(center) = KA * depot - CL*C2
    cmt(eff);
})
print(ode.1c.ka)
```

```R
#> rxode2 2.0.11 model named rx_4caaa6b18411f9babd3e3aafb7840fd4 model (✔ ready). 
#> $state: depot, center
#> $stateExtra: eff
#> $params: V, KA, CL
#> $lhs: C2
```

但是在所有微分方程之前定义的房室是不支持的;所以下面的模型：

```R
ode.1c.ka <- rxode2({
    cmt(eff);
    C2 = center/V;
    d / dt(depot) = -KA * depot
    d/dt(center) = KA * depot - CL*C2
})
```

会出错：

```R
Error in rxModelVars_(obj) : 
  Evaluation error: Compartment 'eff' needs differential equations defined.
```

