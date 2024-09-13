# 使用ubuntu镜像作为基础镜像
FROM ubuntu:20.04
# FROM debian

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y python3 python3-pip && \
    rm -rf /var/lib/apt/lists/*

# 设置工作目录
WORKDIR /app

# 将你的代码和requirements.txt复制到工作目录
COPY . /app

# 安装Python依赖
RUN pip install --no-cache-dir -r requirements.txt

# 运行你的Python脚本
CMD ["python3", "-u", "main.py"]