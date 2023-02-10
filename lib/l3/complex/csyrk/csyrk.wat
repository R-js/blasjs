(module
 (import "env" "matrices" (memory $0 0))
 (export "matrices" (memory $0))
 (export "csyrk" (func $module/csyrk))
 (func $module/csyrk (param $0 i32) (param $1 i32) (param $2 i32) (param $3 i32) (param $4 f64) (param $5 f64) (param $6 f64) (param $7 f64) (param $8 i32)
  (local $9 f64)
  (local $10 i32)
  (local $11 i32)
  (local $12 i32)
  (local $13 f64)
  (local $14 i32)
  (local $15 f64)
  (local $16 i32)
  (local $17 i32)
  (local $18 i32)
  (local $19 i32)
  (local $20 i32)
  (local $21 i32)
  (local $22 f64)
  (local $23 i32)
  (local $24 i32)
  (local $25 i32)
  (local $26 f64)
  (local $27 f64)
  (local $28 i32)
  (local $29 f64)
  (local $30 f64)
  local.get $7
  f64.const 0
  f64.eq
  local.tee $11
  local.get $6
  f64.const 0
  f64.eq
  i32.and
  local.set $19
  i32.const 0
  local.get $6
  f64.const 1
  f64.eq
  local.get $11
  i32.and
  local.tee $21
  i32.const 0
  local.get $3
  local.get $5
  f64.const 0
  f64.eq
  local.get $4
  f64.const 0
  f64.eq
  i32.and
  local.tee $11
  select
  select
  i32.const 1
  local.get $2
  select
  if
   return
  end
  local.get $2
  i32.const 1
  i32.shl
  local.set $16
  local.get $11
  if
   loop $for-loop|0
    local.get $2
    local.get $12
    i32.gt_u
    if
     local.get $10
     local.get $12
     i32.const 1
     i32.add
     i32.const 1
     i32.shl
     i32.add
     local.get $10
     local.get $16
     i32.add
     local.get $0
     select
     i32.const 3
     i32.shl
     local.tee $3
     local.get $10
     local.get $10
     local.get $12
     i32.const 1
     i32.shl
     i32.add
     local.get $0
     select
     i32.const 3
     i32.shl
     local.tee $1
     i32.sub
     local.set $11
     local.get $19
     if
      local.get $8
      if
       local.get $14
       i32.const 0
       local.get $11
       memory.fill $0
       local.get $11
       local.get $14
       i32.add
       local.set $14
      else
       local.get $1
       i32.const 0
       local.get $11
       memory.fill $0
      end
     else
      loop $for-loop|1
       local.get $1
       local.get $3
       i32.lt_u
       if
        local.get $14
        local.get $1
        local.get $8
        select
        local.tee $11
        f64.load $0
        local.set $4
        local.get $11
        local.get $6
        local.get $4
        f64.mul
        local.get $7
        local.get $11
        i32.const 8
        i32.add
        f64.load $0
        local.tee $5
        f64.mul
        f64.sub
        f64.store $0
        local.get $11
        local.get $6
        local.get $5
        f64.mul
        local.get $7
        local.get $4
        f64.mul
        f64.add
        f64.store $0 offset=8
        local.get $1
        i32.const 16
        i32.add
        local.set $1
        local.get $14
        i32.const 16
        i32.add
        local.set $14
        br $for-loop|1
       end
      end
     end
     local.get $12
     i32.const 1
     i32.add
     local.set $12
     local.get $10
     local.get $16
     i32.add
     local.set $10
     br $for-loop|0
    end
   end
   return
  end
  local.get $2
  local.get $2
  i32.const 1
  i32.add
  i32.mul
  local.get $2
  local.get $2
  i32.mul
  i32.const 1
  i32.shl
  local.get $8
  select
  i32.const 3
  i32.shl
  local.set $20
  local.get $1
  if
   local.get $3
   i32.const 1
   i32.shl
   local.set $23
   loop $for-loop|5
    local.get $2
    local.get $12
    i32.gt_u
    if
     local.get $10
     local.get $12
     i32.const 1
     i32.add
     i32.const 1
     i32.shl
     i32.add
     local.get $10
     local.get $16
     i32.add
     local.get $0
     select
     i32.const 3
     i32.shl
     local.set $24
     local.get $10
     local.get $10
     local.get $12
     i32.const 1
     i32.shl
     i32.add
     local.get $0
     select
     i32.const 3
     i32.shl
     local.tee $1
     local.get $16
     i32.const 3
     i32.shl
     i32.rem_u
     local.get $3
     i32.mul
     local.set $11
     loop $for-loop|6
      local.get $1
      local.get $24
      i32.lt_u
      if
       f64.const 0
       local.set $9
       f64.const 0
       local.set $15
       local.get $14
       local.get $1
       local.get $8
       select
       local.tee $25
       f64.load $0
       local.set $26
       local.get $25
       f64.load $0 offset=8
       local.set $27
       i32.const 0
       local.set $17
       loop $for-loop|7
        local.get $17
        local.get $23
        i32.lt_u
        if
         local.get $17
         i32.const 3
         i32.shl
         local.tee $28
         local.get $11
         local.get $20
         i32.add
         i32.add
         local.set $21
         local.get $9
         local.get $18
         local.get $20
         i32.add
         local.get $28
         i32.add
         local.tee $28
         f64.load $0
         local.tee $13
         local.get $21
         f64.load $0
         local.tee $22
         f64.mul
         local.get $28
         f64.load $0 offset=8
         local.tee $29
         local.get $21
         f64.load $0 offset=8
         local.tee $30
         f64.mul
         f64.sub
         f64.add
         local.set $9
         local.get $15
         local.get $13
         local.get $30
         f64.mul
         local.get $29
         local.get $22
         f64.mul
         f64.add
         f64.add
         local.set $15
         local.get $17
         i32.const 2
         i32.add
         local.set $17
         br $for-loop|7
        end
       end
       local.get $4
       local.get $9
       f64.mul
       local.get $5
       local.get $15
       f64.mul
       f64.sub
       local.set $13
       local.get $4
       local.get $15
       f64.mul
       local.get $5
       local.get $9
       f64.mul
       f64.add
       local.set $9
       local.get $19
       i32.eqz
       if
        local.get $13
        local.get $6
        local.get $26
        f64.mul
        local.get $7
        local.get $27
        f64.mul
        f64.sub
        f64.add
        local.set $13
        local.get $9
        local.get $6
        local.get $27
        f64.mul
        local.get $7
        local.get $26
        f64.mul
        f64.add
        f64.add
        local.set $9
       end
       local.get $25
       local.get $13
       f64.store $0
       local.get $25
       local.get $9
       f64.store $0 offset=8
       local.get $1
       i32.const 16
       i32.add
       local.set $1
       local.get $11
       local.get $23
       i32.const 3
       i32.shl
       i32.add
       local.set $11
       local.get $14
       i32.const 16
       i32.add
       local.set $14
       br $for-loop|6
      end
     end
     local.get $12
     i32.const 1
     i32.add
     local.set $12
     local.get $10
     local.get $16
     i32.add
     local.set $10
     local.get $18
     local.get $23
     i32.const 3
     i32.shl
     i32.add
     local.set $18
     br $for-loop|5
    end
   end
  else
   loop $for-loop|2
    local.get $2
    local.get $10
    i32.gt_u
    if
     local.get $12
     local.get $10
     i32.const 1
     i32.add
     i32.const 1
     i32.shl
     i32.add
     local.get $12
     local.get $16
     i32.add
     local.get $0
     select
     i32.const 3
     i32.shl
     local.set $23
     local.get $12
     local.get $12
     local.get $10
     i32.const 1
     i32.shl
     i32.add
     local.get $0
     select
     i32.const 3
     i32.shl
     local.set $1
     loop $for-loop|3
      local.get $1
      local.get $23
      i32.lt_u
      if
       local.get $1
       local.get $16
       i32.const 3
       i32.shl
       i32.rem_u
       local.set $11
       local.get $10
       i32.const 4
       i32.shl
       local.set $17
       f64.const 0
       local.set $9
       f64.const 0
       local.set $15
       i32.const 0
       local.set $18
       loop $for-loop|4
        local.get $3
        local.get $18
        i32.gt_u
        if
         local.get $9
         local.get $11
         local.get $20
         i32.add
         local.tee $24
         f64.load $0
         local.tee $13
         local.get $17
         local.get $20
         i32.add
         local.tee $25
         f64.load $0
         local.tee $22
         f64.mul
         local.get $24
         f64.load $0 offset=8
         local.tee $26
         local.get $25
         f64.load $0 offset=8
         local.tee $27
         f64.mul
         f64.sub
         f64.add
         local.set $9
         local.get $15
         local.get $13
         local.get $27
         f64.mul
         local.get $26
         local.get $22
         f64.mul
         f64.add
         f64.add
         local.set $15
         local.get $18
         i32.const 1
         i32.add
         local.set $18
         local.get $16
         i32.const 3
         i32.shl
         local.tee $24
         local.get $11
         i32.add
         local.set $11
         local.get $17
         local.get $24
         i32.add
         local.set $17
         br $for-loop|4
        end
       end
       local.get $4
       local.get $9
       f64.mul
       local.get $5
       local.get $15
       f64.mul
       f64.sub
       local.set $22
       local.get $4
       local.get $15
       f64.mul
       local.get $5
       local.get $9
       f64.mul
       f64.add
       local.set $9
       local.get $14
       local.get $1
       local.get $8
       select
       local.tee $11
       local.get $19
       if (result f64)
        local.get $22
       else
        local.get $11
        f64.load $0
        local.set $13
        local.get $11
        f64.load $0 offset=8
        local.set $15
        local.get $9
        local.get $21
        if (result f64)
         local.get $15
        else
         local.get $6
         local.get $15
         f64.mul
         local.get $7
         local.get $13
         f64.mul
         f64.add
         local.set $9
         local.get $6
         local.get $13
         f64.mul
         local.get $7
         local.get $15
         f64.mul
         f64.sub
         local.set $13
         local.get $9
        end
        f64.add
        local.set $9
        local.get $22
        local.get $13
        f64.add
       end
       f64.store $0
       local.get $11
       local.get $9
       f64.store $0 offset=8
       local.get $1
       i32.const 16
       i32.add
       local.set $1
       local.get $14
       i32.const 16
       i32.add
       local.set $14
       br $for-loop|3
      end
     end
     local.get $10
     i32.const 1
     i32.add
     local.set $10
     local.get $12
     local.get $16
     i32.add
     local.set $12
     br $for-loop|2
    end
   end
  end
 )
)