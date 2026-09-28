# 只读状态网页

只依赖 Python 3.10+ 标准库，不需要安装 Node、数据库或其他包。不启动或重启实验，读取现有 `results/full-paper`。Linux/WSL 上显示 CPU、内存、交换内存、磁盘；每 5 秒刷新。实验详情包含最新已写完整的时间步、Newton 迭代数、det F 和最近日志。缓存/写缓冲可能使时间步显示略有滞后。

## 在现有服务器添加网页

先检查 `git status`。若服务器 Codex 已修改文件，请保留这些修改并合并更新，禁止 reset --hard 或删除结果。干净工作区可以：

```sh
cd /home/karen/src/dealii-research
git pull --ff-only
python3 tests/test_dashboard.py
mkdir -p results/dashboard
nohup python3 scripts/dashboard.py --host 0.0.0.0 --port 8766 \
  --token-file "$HOME/.config/dealii-dashboard/token" \
  > results/dashboard/server.log 2>&1 < /dev/null &
```

先检查 8766 是否已有服务，已有该监控进程就复用，不重复启动。token 文件第一次自动创建，权限 600；不要提交到仓库或转发给其他人。查看 token 后，在服务器浏览器打开 `http://localhost:8766/?token=你的token`。日志不会记录 token。Linux 本地请求可验证：

```sh
curl -fsS -H "Authorization: Bearer $(cat "$HOME/.config/dealii-dashboard/token")" \
  http://127.0.0.1:8766/api/status
```

Windows 到 WSL 的 localhost 转发若可用，也可在 Windows 浏览器使用该地址。

## 从 Mac 访问

推荐通过现有 SSH 加密转发，无需增加 Windows 防火墙开放端口：

```sh
ssh -N -L 18766:127.0.0.1:8766 -b 10.6.214.71 \
  -o IPQoS=none -o KexAlgorithms=curve25519-sha256 karen@10.6.214.78
```

然后打开 `http://localhost:18766/?token=你的token`，保持该连接。若网络地址变化，更新相应 IP。

若希望 Mac 直接打开 `http://10.6.214.78:8766/?token=你的token`，需由服务器端将 Windows TCP 8766 转发到**当前 WSL IP**的 TCP 8766，并将防火墙访问范围限定为你的 Mac 地址。已有的 22 端口转发不会自动转发 8766。直接 HTTP 没有传输加密，优先使用 SSH 转发；不要开放到公网。

## 状态判定与边界

- completed 必须来自 `summary.json` 的 completed=true；这不代表科学结论已验证。
- running 必须发现对应输出路径的 giesekus 进程；只存在目录不会视作运行。
- failed 来自算例 status.json；目录已有内容但无运行进程、无完成记录时显示“未完成 / 待检查”。
- 主流程按 pipeline.pid 与进程命令核对。其他恢复脚本接管时可能不显示主流程存活，应结合阶段、算例进程和日志判断。
- 阶段显示的是 phase.txt 原文和写入时间；连接失败时明确提示当前页面为旧数据。
- CPU 是全机逻辑 CPU 的总体利用率，内存是全机用量，不能归因于这一个项目。GPU 不使用。
- 网页是观察工具，不会主动修复、发送通知或自动重启进程。停止网页不影响计算，计算停止也不会停止网页。
- 通过 `--root /absolute/path/to/results` 可观察其他结果根目录；不会扫描用户家目录。日志只允许固定名称及清单内算例，最多读取尾部 24 KB；无文件浏览或命令执行接口。
