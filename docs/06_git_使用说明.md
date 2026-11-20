# Git 使用说明

> 给不熟悉 git 的使用者。**日常只有三条命令**：`git add -A`、`git commit -m "..."`、`git push`。

## 1. 现在已经做好的部分

| 项 | 状态 |
|---|---|
| 本地仓库 | ✅ 已初始化，分支 `main`，2 个提交 |
| 提交身份 | `MartinYang-SCP <MartinYang-SCP@users.noreply.github.com>`（仓库级配置，可改） |
| 远程仓库 | ✅ `origin` = `git@github.com:MartinYang-SCP/ScaffoldAssistedDNASynthesisTool.git` |
| `.gitignore` | ✅ 已忽略 `__pycache__/`、`build/`、测试缓存、本机密钥 |
| SSH 密钥 | ✅ 已生成 `gitkey_ed25519`（私钥）+ `sshkey_for_github.txt`（公钥文本），**都被 ignore** |

**唯一没做的**：把提交推送到 GitHub。原因见下节——不是 git 用错，是这台机器上两条路都不通。

## 2. 为什么没有自动推上去

| 通道 | 现象 | 结论 |
|---|---|---|
| HTTPS (`github.com:443`) | 连续超时（`Failed to connect after 21s`） | 网络层面不通 |
| SSH（`git@github.com`） | 握手能成功，但 git 调 ssh 时要经过它自带的 `sh.exe`，而 `sh.exe` 在当前受限环境里起不来（`couldn't create signal pipe, Win32 error 5`） | 只能在普通终端里执行 |

也就是说：**在你自己打开的 PowerShell / Git Bash 里跑一次 `git push` 就行**，那边两个限制都不存在。

## 3. 一次性设置：把公钥加到 GitHub

1. 打开 `sshkey_for_github.txt`，复制里面**一整行**（`ssh-ed25519 AAAA... scaffold-tool@dsh`）；
2. 打开 <https://github.com/settings/keys> → **New SSH key**；
3. Title 随便写（例如 `work-pc`），Key 粘贴，保存；
4. 在仓库目录里验证：

```powershell
ssh -T git@github.com
# 期望看到：Hi MartinYang-SCP! You've successfully authenticated...
```

> 想更"标准"一点：把 `gitkey_ed25519` 和 `.pub` 移到 `%USERPROFILE%\.ssh\id_ed25519`，
> 然后 `git config --unset core.sshCommand`。不改也能用，仓库里已经配好 `core.sshCommand` 指向它。

## 4. 推送（复制粘贴即可）

在**普通终端**（不是受限环境）里，进入项目目录：

```powershell
cd "D:\THU\大四上\核酸纳米材料\小课题"

# 远端可能已经有 1 个提交（建仓时自动生成的 README），先取回来合并
git fetch origin
if (git rev-parse --verify origin/main 2>$null) {
    git merge --allow-unrelated-histories -X ours origin/main -m "merge remote initial commit"
}

git push -u origin main
```

如果远端是空的，`if` 那段会自动跳过，直接 push。

> `-X ours` 的意思是：万一远端那个自动生成的 README 与本项目的 README 撞名，
> 保留**本地**版本。远端那个提交仍然留在历史里，不会被抹掉。

## 5. 日常改动怎么提交

```powershell
git status                  # 看改了什么
git add -A                  # 全部加入暂存区
git commit -m "说明这次改了什么"
git push                    # 推到 GitHub
```

提交信息写清楚"为什么改"，比写"update"有用得多。

## 6. 常见问题

**Q：`git push` 报 `rejected (fetch first)`**
A：远端有你本地没有的提交。先 `git pull --rebase`，再 `git push`。

**Q：改了文件后 `git status` 看不到**
A：多半被 `.gitignore` 忽略了（`build/`、`__pycache__/` 本来就不该提交）。

**Q：不小心提交了密钥**
A：立刻删文件、`git rm --cached 文件`、提交，并到 GitHub 撤销那把 key。本项目已用 `.gitignore` 屏蔽 `gitkey_*`、`sshkey_for_github.txt`、`.ssh_known_hosts`。

**Q：想撤销上一次提交但保留改动**
A：`git reset --soft HEAD~1`。

**Q：中文文件名显示成八进制转义**
A：已设置 `core.quotepath=false`，正常显示中文。
