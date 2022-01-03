# Adaptive-Speech-Dereverberation-based-on-QR-MCLP-model
根据论文《Multi-Channel Linear Prediction Speech Dereverberation Algorithm Based on QR-RLS Adaptive Filter》写了一个C代码去实现双通道语音去混响的功能。

代码是针对双通道语音32k进行去混响的，一开始先将32k语音进行正交分解，分解成低16k的成分和高16k成分。然后直接对低频部分（双通道16k）语音进行QR-MCLP模型去混响。去混响结束后直接与高频部分进行拼接得到输出。

由于语音以及混响成分主要集中在低频部分，所以高频部分不做任何处理依然可以有比较好的性能效果。

目前算法已经在仿真和实录音频上测试过，效果良好。

