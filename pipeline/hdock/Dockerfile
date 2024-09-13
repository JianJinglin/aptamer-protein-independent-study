# 使用ubuntu镜像作为基础镜像
FROM ubuntu:20.04
# FROM debian

# 设置工作目录
WORKDIR /app

# 将你的代码和requirements.txt复制到工作目录
COPY . /app

RUN apt-get update
RUN apt-get install -y libfftw3-dev
RUN chmod +x /app/hdock
RUN export PATH=$PATH:/app/hdock

# 运行你的Python脚本
CMD ["./hdock"]