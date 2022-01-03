# Adaptive-Speech-Dereverberation-based-on-QR-MCLP-model
根据论文《Multi-Channel Linear Prediction Speech Dereverberation Algorithm Based on QR-RLS Adaptive Filter》写了一个C程序去实现双通道语音去混响的功能。

代码是针对双通道32kHz语音进行去混响的，一开始先将32kHz语音进行正交分解，分解成低16kHz的成分和高16kHz成分。然后直接对低频部分（双通道16kHz）语音进行QR-MCLP模型去混响。去混响结束后直接与高频部分进行拼接得到输出。

由于语音以及混响成分主要集中在低频部分，所以高频部分不做任何处理依然可以有比较好的性能效果。

模型的主体是一个MCLP循环最小二乘自适应滤波器，使用在求解参数时使用QR分解来保证模型的鲁棒性。QR分解以及回溯部分运行时间较长。

目前算法已经在仿真和实录音频上测试过，效果良好。

需要注意的是音频的输入输出都是PCM格式的，如果是WAV格式或其他需要进行转换。

论文部分已附在工程文件中。

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

translate:

According to the paper "Multi-Channel Linear Prediction Speech Dereverberation Algorithm Based on QR-RLS Adaptive Filter", this is a C program to realize the function of dual-channel speech de-reverberation.

The code is for 32kHz dual-channel speech de-reverberation. At the beginning, the raw speech is orthogonally decomposed into a low 16kHz component and a high 16kHz component. Then the program directly uses the QR-MCLP model to dereverberate the audio on the low frequency part (16kHz dual channel audio) speech. After the de-reverberation, the processed low-frequency part audio is directly spliced with the high-frequency part to get the output.

Since the voice and reverberation components are mainly concentrated in the low frequency part, we can still have a superb performance effect without any processing on the high frequency part.

The main body of the model is an MCLP Recursive Least Squares adaptive filter, which uses QR decomposition to ensure the robustness of the model. The QR decomposition and backtracking part take a long time to run.

At present, the algorithm has been tested on simulated and recorded audio, and the effect is good.

It should be noted that the audio input and output are all in PCM format,  and the WAV format or other needs to be converted.

The paper has been attached to this repository.

