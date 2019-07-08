# 无脑~超简单本地Git配合GitHub

> 刘小泽写于18.12.24
>
> 花花想用GitHub，于是我写了这个

### 第一步 简单注册

注册GitHub账号，邮箱是之后常用到的

### 第二步 初步配置

#### 先在本地terminal

本地检查是否有Git，看看`git`就知道电脑上有没有了

在本地新建一个文件夹，比如叫`GEO` ，用来存放R代码

```shell
mkdir ~/Git/GEO
cd ~/Git/GEO
git init # 初始化git
# 然后将本地git与GitHub联系起来
cd ~/.ssh # (注意这是隐藏文件夹，用ls -la才能查看)
ssh-keygen -t rsa -C your@mail.com # 改一下邮箱名就好
# 然后看到.ssh文件夹中存在了id_rsa.pub
cat id_rsa.pub #然后将内容复制下来
```

#### 再去浏览器

登录上GitHub账号=》右上角头像=》`Settings` =》`SSH and GPG keys` =》右上角`New SSH keys` =》将之前复制的粘贴上

#### 再回到本地terminal

输入` ssh -T git@github.com` ，如果出现`Hi xxx! You've successfully autheticated` **就成功啦！HooRay！🥰**

### 第三步 上传代码

#### 先在浏览器

先在GitHub上创建一个自己的Repository，很简单的过程

还是先点右上角头像=》`Your Profile` =》然后会看到`Repositories`的选项=》点击然后再点击绿色的`New`=》填写`Repository name` （比如填GEO）=》其他可以不用填，然后底部`Create repository` **就成功啦！HooRay！**🥰

然后会刷新一下，你会看到一个类似于`https://github.com/XXXX/GEO.git` 的链接，复制下来

#### 再在本地

还是在terminal中

```shell
cd ~/Git/GEO
git remote add origin https://github.com/YOUR_NAME/GEO.git 
# 就是刚才复制的链接，这样就把本地和网络端联系起来了
```

然后可以新建你自己的代码了，比如一个或几个关于GEO的R脚本，从其他地方复制到`~/Git/GEO`中

接着，`git add .` 【表示将当前文件夹中的全部新增/新修改的文件准备好】

然后`git status` 【看看刚才的操作增加了哪些文件，是不是自己想要的；如果不是，也有办法去掉某几个`git reset HEAD <FILE>` (这个不重要现在！)】 

然后`git commit -m "你想写的备注"` 

最后`git push -u origin master` 【`-u`参数只需要第一次输入，以后只需要输入`git push origin master`】

### 最后，请注意

本地的一个文件夹如`GEO`只对应GitHub的一个`Repository`：

如果自己本地有多个文件夹，比如还有shell脚本的文件夹，perl脚本的文件夹，**一定要先在GitHub上新建好对应的`Repository`，然后再按第三步重新走一遍**

这样就确保自己的每个文件夹中的代码都能同步到GitHub做备份

>GitHub的重要性不用多说，可以随时记录你的脚本改动，并且可以及时恢复到任何版本
>
>好啦！以上就是超级简单的Git小教程。希望对你有帮助



---

> 日常整理git问题

##### 2019.7.8

1：总是提示：`Enter passphrase for key '/Users/reedliu1/.ssh/id_rsa':`

- 使用`ssh-add ~/.ssh/id_rsa` 然后添加密码就好了

